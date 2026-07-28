"""Generate a region-connectivity OWL file for one connectomics dataset.

For each neuron in a dataset that has a VFB individual, this emits the neuron's
per-neuropil-region synaptic tallies as axiom-annotated relations, matching the
representation already in VFB's production graph (pdb):

  * has presynaptic terminals in  (RO:0002113)  -> annotations: Tbars, downstream
  * has postsynaptic terminal in  (RO:0002110)  -> annotation:  upstream

Counts are pulled fresh from neuprint (neuprint-python) for the neuprint dataset
that matches the VFB Site node. The count annotations use the n2o `custom:`
namespace so that `owl2neo4jcsv` turns each axiom into a Neo4j edge carrying the
counts as edge properties (side-load into pdb via vfb-pipeline-dumps).

Each edge is an ObjectPropertyAssertion to a VFB region *individual*; the
individual is side-specific, so laterality is carried by the target.

Usage:
    python3 make_region_connectivity.py <dataset_key> <neuprint_token_file> <out_dir> <report_dir>
"""
import os
import re
import sys
from urllib.parse import unquote

import pandas as pd
import yaml
import neuprint
from neuprint import NeuronCriteria as NC, fetch_neurons

from cypher_query import query_neo4j, auth

# VFB Neo4j endpoints. pdb (production) holds released datasets + region
# individuals; kb (knowledge base, source) holds datasets not yet in pdb.
ENDPOINTS = {
    "pdb": "https://pdb.virtualflybrain.org/db/neo4j/tx/commit",
    "kb": "http://kb.virtualflybrain.org/db/neo4j/tx/commit",
}
DEFAULT_ENDPOINT = "kb"
# Default for standalone calls; main() overrides per-dataset (see DATASETS).
KB = ENDPOINTS[DEFAULT_ENDPOINT]

# --- Relations & count annotations -----------------------------------------
RO_POST = "http://purl.obolibrary.org/obo/RO_0002110"  # has postsynaptic terminal in
RO_PRE = "http://purl.obolibrary.org/obo/RO_0002113"   # has presynaptic terminals in
CUSTOM = "http://n2o.neo/custom/"                       # -> Neo4j edge properties

# neuprint roi_counts_df column -> custom annotation local name, per relation.
# (matches existing pdb edge properties)
# Each relation has a PRIMARY count that gates edge emission — the edge means the
# neuron has those *terminals* in the region, so a presynaptic edge requires
# T-bars (pre), not just a stray downstream connection count — plus optional
# secondary counts. (relation, primary neuprint col, primary local,
# ((secondary col, secondary local), ...))
EDGE_SPECS = (
    (RO_POST, "upstream", "upstream", ()),
    (RO_PRE, "pre", "Tbars", (("downstream", "downstream"),)),
)
# deterministic ordering of count annotations within an axiom
_LOCAL_ORDER = {"upstream": 0, "Tbars": 1, "downstream": 2}

# Each edge carries only the merged total per count. Region individuals are
# side-specific, so laterality is recoverable downstream from the region's side
# (via part_of in the kb) joined with the neuron's soma side (has_soma_location).

# Datatype of the count annotations. pdb stores these as float
# (neo4j "LIST OF FLOAT", e.g. [545.0]); we match that. Reports cast the counts
# back to integer for readability.
COUNT_DATATYPE = "http://www.w3.org/2001/XMLSchema#float"

# number of axioms written to the per-file QC sample (shape-checked at build time,
# since the full file is too large for robot verify)
SAMPLE_AXIOMS = 5000

# --- Per-dataset configuration ---------------------------------------------
# Datasets are configured in src/ontology/dataset_config.yaml (shared with the
# Makefile); see that file for the meaning of each key.
DATASET_CONFIG = os.path.join(
    os.path.dirname(__file__), "..", "ontology", "dataset_config.yaml"
)


def load_datasets(path=DATASET_CONFIG):
    with open(path) as fh:
        return yaml.safe_load(fh)


DATASETS = load_datasets()


def get_neuprint_dataset(site, url=KB):
    """Derive the neuprint dataset string (e.g. 'hemibrain:v1.2.1') from the
    Site node, so the counts come from the version VFB is pinned to."""
    df = query_neo4j(
        'MATCH (s:Site {short_form:"%s"}) RETURN s.link_base AS lb, s.homepage AS hp'
        % site,
        url=url,
    )
    if df is False or df.empty:
        raise RuntimeError(f"Site node {site} not found")
    for col in ("lb", "hp"):
        val = df.iloc[0][col]
        text = val[0] if isinstance(val, list) and val else val
        if isinstance(text, str):
            m = re.search(r"dataset=([^&\s]+)", text)
            if m:
                return unquote(m.group(1))  # link_base may be URL-encoded (%3A -> :)
    raise RuntimeError(f"Could not derive neuprint dataset from Site {site}")


def get_bodyid_to_neuron(site, url=KB):
    """bodyId -> VFB neuron individual IRI (non-deprecated).

    Uses the :Individual label (present in both pdb and the kb), since the
    derived :Neuron label only exists in pdb.
    """
    df = query_neo4j(
        'MATCH (n:Individual)-[r:database_cross_reference|hasDbXref]->'
        '(s:Site {short_form:"%s"}) '
        "WHERE (NOT EXISTS(n.deprecated)) OR (NOT n.deprecated[0]) "
        "RETURN toInteger(r.accession[0]) AS bodyId, n.iri AS iri" % site,
        url=url,
    )
    df = df.dropna(subset=["bodyId"])
    df["bodyId"] = df["bodyId"].astype("int64")
    return dict(zip(df["bodyId"], df["iri"]))


def get_roi_to_region(roi_dataset, url=KB):
    """neuprint ROI name -> VFB region-individual IRI.

    The ROI name is the region label with the ' on <template>' suffix stripped,
    which equals the neuprint ROI name.
    """
    df = query_neo4j(
        "MATCH (reg)-[:has_source]->(:DataSet {short_form:\"%s\"}) "
        "RETURN reg.iri AS iri, reg.label AS label" % roi_dataset,
        url=url,
    )
    return {_roi_from_label(row["label"]): row["iri"] for _, row in df.iterrows()}


def _roi_from_label(label):
    """VFB region label -> neuprint ROI name.

    Strips the ' on <template>' suffix, and removes literal backslashes: some
    VFB labels store an escaped apostrophe (e.g. "a\\'L(R)") while neuprint ROI
    names use a plain apostrophe, so MB alpha'/beta' lobe ROIs would otherwise
    fail to map and be dropped.
    """
    if isinstance(label, list):
        label = label[0]
    return label.split(" on ")[0].replace("\\", "")


def fetch_roi_counts(client, body_ids):
    """neuprint roi_counts_df restricted to the bodyIds we have VFB IRIs for."""
    frames = []
    ids = list(body_ids)
    chunk = 5000
    for i in range(0, len(ids), chunk):
        _, roidf = fetch_neurons(NC(bodyId=ids[i : i + chunk]), client=client)
        if roidf is not None and not roidf.empty:
            frames.append(roidf)
    if not frames:
        return pd.DataFrame(columns=["bodyId", "roi", "pre", "post", "downstream", "upstream"])
    return pd.concat(frames, ignore_index=True)


def roi_groups_from_hierarchy(root):
    """neuprint ROI name -> region group ('brain' or 'vnc') from the ROI
    hierarchy: any 'VNC' node (and its whole subtree) is 'vnc', the rest 'brain'.

    Walks from the root, so it handles a VNC-only dataset whose root *is* 'VNC'
    (manc), a whole-CNS dataset with a 'VNC' branch (male-cns), and a brain-only
    dataset with no 'VNC' node (hemibrain -> all 'brain')."""
    groups = {}

    def walk(node, g):
        if node["name"] == "VNC":
            g = "vnc"
        groups[node["name"]] = g
        for ch in node.get("children") or []:
            walk(ch, g)

    walk(root, "brain")
    return groups


def get_roi_groups(client):
    return roi_groups_from_hierarchy(client.meta["roiHierarchy"])


# --- OWL (functional syntax) emission --------------------------------------
def _lit(n):
    return '"%d"^^<%s>' % (int(n), COUNT_DATATYPE)


# trailing (R)/(L) or (R<n>)/(L<n>) — the latter for protocerebral-bridge halves
def _edge_counts(row, primary_col, primary_local, secondary):
    """Counts for one relation on a neuprint row, or None if the primary
    (terminal) count is absent/zero — so no edge is emitted from a secondary
    count alone. No synapse-count threshold beyond >= 1."""
    pv = row.get(primary_col)
    if not (pd.notna(pv) and pv > 0):
        return None
    counts = {primary_local: pv}
    for col, local in secondary:
        v = row.get(col)
        if pd.notna(v) and v > 0:
            counts[local] = v
    return counts


def _count_annotations(counts):
    return [
        "Annotation(<%s%s> %s)" % (CUSTOM, local, _lit(counts[local]))
        for local in sorted(counts, key=_LOCAL_ORDER.get)
    ]


def _relation_axiom(neuron_iri, target_iri, prop, anns):
    return "ObjectPropertyAssertion(%s <%s> <%s> <%s>)" % (
        " ".join(anns), prop, neuron_iri, target_iri,
    )


def build_axioms(counts, bodyid2neuron, roi2region):
    """One ObjectPropertyAssertion per (neuron, region individual, relation).
    Region individuals are side-specific, so laterality is carried by the target."""
    axioms = []
    used_props = set()
    # drop rows for ROIs with no VFB region up front (e.g. male-cns's ~5k
    # optic-lobe column ROIs) so we don't iterate millions of skipped rows
    counts = counts[counts["roi"].isin(roi2region)]
    for _, row in counts.iterrows():
        neuron = bodyid2neuron.get(int(row["bodyId"]))
        target = roi2region.get(row["roi"])
        if neuron is None or target is None:
            continue
        for rel, pcol, plocal, secondary in EDGE_SPECS:
            cs = _edge_counts(row, pcol, plocal, secondary)
            if cs is not None:
                axioms.append(_relation_axiom(neuron, target, rel, _count_annotations(cs)))
                used_props.add(rel)
    return axioms, used_props


def write_ontology(path, dataset_key, neuprint_dataset, axioms, used_props):
    onto_iri = (
        "http://virtualflybrain.org/data/VFB/OWL/region_connectivity_%s.owl"
        % dataset_key
    )
    # declare exactly the custom annotation properties actually used
    used_anns = set(re.findall(r"<http://n2o\.neo/custom/([A-Za-z_]+)>", "\n".join(axioms)))
    decls = [
        "Declaration(AnnotationProperty(<%s%s>))" % (CUSTOM, n)
        for n in sorted(used_anns)
    ]
    decls += ["Declaration(ObjectProperty(<%s>))" % p for p in sorted(used_props)]
    # neuprint dataset is "<name>:<version>", e.g. hemibrain:v1.2.1
    version = neuprint_dataset.partition(":")[2] or neuprint_dataset
    lines = [
        "Prefix(:=<http://purl.obolibrary.org/obo/>)",
        "Prefix(owl:=<http://www.w3.org/2002/07/owl#>)",
        "Prefix(xsd:=<http://www.w3.org/2001/XMLSchema#>)",
        "Prefix(dcterms:=<http://purl.org/dc/terms/>)",
        "Prefix(custom:=<%s>)" % CUSTOM,
        "",
        "Ontology(<%s>" % onto_iri,
        'Annotation(rdfs:comment "Region (neuropil) synaptic connectivity for '
        '%s, generated from neuprint %s.")' % (dataset_key, neuprint_dataset),
        'Annotation(dcterms:title "VFB region connectivity: %s")' % dataset_key,
        'Annotation(dcterms:identifier "%s")' % dataset_key,
        'Annotation(dcterms:source "%s")' % neuprint_dataset,
        'Annotation(owl:versionInfo "%s")' % version,
        "",
    ]
    lines += decls + [""] + axioms + [")"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_roi_report(path, counts, roi2region, roi_group=None):
    """Report every neuprint ROI seen in the data: its VFB region individual (or
    blank if unmapped), its `group`, and how many neurons touch it. The group is
    `unmapped` for any ROI with no region individual (these are dropped from the
    OWL — typically neuprint super-ROIs like SNP(R)/INP or optic-lobe columns; a
    large unmapped ROI would signal a missing region individual), otherwise the
    region's brain/vnc group."""
    roi_group = roi_group or {}
    per_roi = counts.groupby("roi")["bodyId"].nunique()
    rows = [(roi,
             "unmapped" if roi not in roi2region else roi_group.get(roi, ""),
             roi2region.get(roi, ""),
             int(neurons))
            for roi, neurons in per_roi.items()]
    # unmapped ROIs first (they need attention), then alphabetical
    rows.sort(key=lambda r: (r[2] != "", r[0]))
    with open(path, "w") as fh:
        fh.write("roi\tgroup\ttarget\tneurons\n")
        for roi, grp, iri, neurons in rows:
            fh.write("%s\t%s\t%s\t%d\n" % (roi, grp, iri, neurons))
    unmapped = [r for r in rows if r[2] == ""]
    print("ROI mapping: %d ROIs, %d unmapped -> %s" % (len(rows), len(unmapped), path))
    if unmapped:
        print("  unmapped:", ", ".join(r[0] for r in unmapped))
    return rows


def main():
    # args: <dataset> <token_file> <out_dir> <report_dir>
    dataset_key, token_file, out_dir, report_dir = sys.argv[1:5]
    if dataset_key not in DATASETS:
        sys.exit("unknown dataset %r; configured: %s"
                 % (dataset_key, ", ".join(sorted(DATASETS))))
    cfg = DATASETS[dataset_key]
    token = open(os.path.expanduser(token_file)).read().strip()

    # Default to the kb: dataset updates land there first (needed to migrate to
    # a new dataset version before it reaches pdb).
    url = ENDPOINTS[cfg.get("endpoint", "kb")]

    neuprint_dataset = get_neuprint_dataset(cfg["site"], url)
    print("neuprint dataset:", neuprint_dataset)

    bodyid2neuron = get_bodyid_to_neuron(cfg["site"], url)
    print("VFB neurons:", len(bodyid2neuron))

    roi2region = get_roi_to_region(cfg["roi_dataset"], url)
    print("regions mapped:", len(roi2region))

    client = neuprint.Client("https://neuprint.janelia.org", dataset=neuprint_dataset, token=token)
    roi_group = get_roi_groups(client)
    counts = fetch_roi_counts(client, bodyid2neuron.keys())
    print("neuprint roi-count rows:", len(counts))

    os.makedirs(report_dir, exist_ok=True)
    # roi_group only annotates the report's `group` column (brain/vnc); one file
    # per dataset (a whole-CNS dataset stays well under GitHub's size limit).
    write_roi_report(os.path.join(report_dir, "roi_mapping.tsv"), counts, roi2region, roi_group)

    axioms, used_props = build_axioms(counts, bodyid2neuron, roi2region)
    print("axioms:", len(axioms))

    qc_dir = os.path.join(out_dir, "qc")
    os.makedirs(qc_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "region_connectivity_%s.ofn" % dataset_key)
    write_ontology(out_path, dataset_key, neuprint_dataset, axioms, used_props)
    # small sample for build-time QC (full file is too large for robot verify)
    write_ontology(os.path.join(qc_dir, "region_connectivity_%s.ofn" % dataset_key),
                   dataset_key, neuprint_dataset, axioms[:SAMPLE_AXIOMS], used_props)
    print("wrote", out_path)


if __name__ == "__main__":
    main()

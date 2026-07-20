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

A region is emitted as either:
  * an ObjectPropertyAssertion to a VFB region *individual* (datasets that have
    per-dataset painted-domain individuals, e.g. hemibrain), or
  * a ClassAssertion over ObjectSomeValuesFrom to an FBbt *class* (datasets with
    no region individuals).

Usage:
    python3 make_region_connectivity.py <dataset_key> <neuprint_token_file> <out.owl.ofn>
"""
import os
import re
import sys

import pandas as pd
import neuprint
from neuprint import NeuronCriteria as NC, fetch_neurons

from cypher_query import query_neo4j

# pdb is the read endpoint holding the VFB individuals + region individuals.
PDB = "https://pdb.virtualflybrain.org/db/neo4j/tx/commit"

# --- Relations & count annotations -----------------------------------------
RO_POST = "http://purl.obolibrary.org/obo/RO_0002110"  # has postsynaptic terminal in
RO_PRE = "http://purl.obolibrary.org/obo/RO_0002113"   # has presynaptic terminals in
CUSTOM = "http://n2o.neo/custom/"                       # -> Neo4j edge properties

# neuprint roi_counts_df column -> custom annotation local name, per relation.
# (matches existing pdb edge properties)
POST_ANNOTATIONS = {"upstream": "upstream"}
PRE_ANNOTATIONS = {"pre": "Tbars", "downstream": "downstream"}

# Datatype of the count annotations. pdb stores these as float
# (neo4j "LIST OF FLOAT", e.g. [545.0]); we match that. Reports cast the counts
# back to integer for readability.
COUNT_DATATYPE = "http://www.w3.org/2001/XMLSchema#float"

# --- Per-dataset configuration ---------------------------------------------
# site        : VFB Site node short_form (the connectome data source).
# roi_dataset : VFB DataSet short_form whose members are the region *individuals*
#               (painted domains). Set to None for datasets whose regions must map
#               to FBbt *classes* instead.
DATASETS = {
    "hemibrain": {
        "site": "neuprint_JRC_Hemibrain_1point2point1",
        "roi_dataset": "Xu2020roi",
    },
}


def get_neuprint_dataset(site):
    """Derive the neuprint dataset string (e.g. 'hemibrain:v1.2.1') from the
    Site node, so the counts come from the version VFB is pinned to."""
    df = query_neo4j(
        'MATCH (s:Site {short_form:"%s"}) RETURN s.link_base AS lb, s.homepage AS hp'
        % site,
        url=PDB,
    )
    if df is False or df.empty:
        raise RuntimeError(f"Site node {site} not found in pdb")
    for col in ("lb", "hp"):
        val = df.iloc[0][col]
        text = val[0] if isinstance(val, list) and val else val
        if isinstance(text, str):
            m = re.search(r"dataset=([^&\s]+)", text)
            if m:
                return m.group(1)
    raise RuntimeError(f"Could not derive neuprint dataset from Site {site}")


def get_bodyid_to_neuron(site):
    """bodyId -> VFB neuron individual IRI (non-deprecated)."""
    df = query_neo4j(
        'MATCH (n:Neuron)-[r:database_cross_reference|hasDbXref]->'
        '(s:Site {short_form:"%s"}) '
        "WHERE (NOT EXISTS(n.deprecated)) OR (NOT n.deprecated[0]) "
        "RETURN toInteger(r.accession[0]) AS bodyId, n.iri AS iri" % site,
        url=PDB,
    )
    df = df.dropna(subset=["bodyId"])
    df["bodyId"] = df["bodyId"].astype("int64")
    return dict(zip(df["bodyId"], df["iri"]))


def get_roi_to_region(roi_dataset):
    """ROI name -> (target IRI, kind) using VFB region individuals.

    kind is 'individual'. The ROI name is the region label with the
    ' on <template>' suffix stripped, which equals the neuprint ROI name.
    """
    df = query_neo4j(
        "MATCH (reg)-[:has_source]->(:DataSet {short_form:\"%s\"}) "
        "RETURN reg.iri AS iri, reg.label AS label" % roi_dataset,
        url=PDB,
    )
    mapping = {}
    for _, row in df.iterrows():
        mapping[_roi_from_label(row["label"])] = (row["iri"], "individual")
    return mapping


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


def load_roi_mapping_tsv(path):
    """ROI name -> (target IRI, kind) from a curated TSV for datasets whose
    regions map to FBbt classes (or a mix). Columns: roi, iri, kind.

    A good source for building such a TSV is EM_neuropil_synonyms.owl in the
    drosophila-anatomy-developmental-ontology repo, which maps FBbt neuropil
    classes to ROI names via oboInOwl:hasRelatedSynonym (e.g. FBbt_00045032 <-
    "SIP(R)"/"SIP_R")."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return {
        r["roi"]: (_obo_iri(r["iri"]), r.get("kind") or "class")
        for _, r in df.iterrows()
    }


def _obo_iri(curie_or_iri):
    if curie_or_iri.startswith("http"):
        return curie_or_iri
    prefix, local = curie_or_iri.split(":", 1)
    return f"http://purl.obolibrary.org/obo/{prefix}_{local}"


def fetch_roi_counts(neuprint_dataset, token, body_ids):
    """neuprint roi_counts_df restricted to the bodyIds we have VFB IRIs for."""
    client = neuprint.Client(
        "https://neuprint.janelia.org", dataset=neuprint_dataset, token=token
    )
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


# --- OWL (functional syntax) emission --------------------------------------
def _lit(n):
    return '"%d"^^<%s>' % (int(n), COUNT_DATATYPE)


def _annotations(row, spec):
    # No synapse-count threshold: any positive count (>= 1) is kept, matching the
    # existing VFB data. To filter, raise the `> 0` gate to a `min_count`.
    return [
        "Annotation(<%s%s> %s)" % (CUSTOM, local, _lit(row[col]))
        for col, local in spec.items()
        if pd.notna(row.get(col)) and row[col] > 0
    ]


def _relation_axiom(neuron_iri, target_iri, kind, prop, anns):
    ann = " ".join(anns)
    if kind == "individual":
        return "ObjectPropertyAssertion(%s <%s> <%s> <%s>)" % (
            ann, prop, neuron_iri, target_iri,
        )
    # class target: existential class assertion
    return "ClassAssertion(%s ObjectSomeValuesFrom(<%s> <%s>) <%s>)" % (
        ann, prop, target_iri, neuron_iri,
    )


def build_axioms(counts, bodyid2neuron, roi2region):
    axioms = []
    used_props = set()
    for _, row in counts.iterrows():
        neuron = bodyid2neuron.get(int(row["bodyId"]))
        region = roi2region.get(row["roi"])
        if neuron is None or region is None:
            continue
        target_iri, kind = region
        # has postsynaptic terminal in (inputs)
        post_anns = _annotations(row, POST_ANNOTATIONS)
        if post_anns:
            axioms.append(_relation_axiom(neuron, target_iri, kind, RO_POST, post_anns))
            used_props.add(RO_POST)
        # has presynaptic terminals in (outputs)
        pre_anns = _annotations(row, PRE_ANNOTATIONS)
        if pre_anns:
            axioms.append(_relation_axiom(neuron, target_iri, kind, RO_PRE, pre_anns))
            used_props.add(RO_PRE)
    return axioms, used_props


def write_ontology(path, dataset_key, neuprint_dataset, axioms, used_props):
    onto_iri = (
        "http://virtualflybrain.org/data/VFB/OWL/region_connectivity_%s.owl"
        % dataset_key
    )
    decls = [
        "Declaration(AnnotationProperty(<%s%s>))" % (CUSTOM, n)
        for n in sorted(set(POST_ANNOTATIONS.values()) | set(PRE_ANNOTATIONS.values()))
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


def write_roi_report(path, counts, roi2region):
    """Report every neuprint ROI seen in the data, what it maps to (a VFB region
    individual or an FBbt class) and how many neurons touch it. ROIs with no VFB
    mapping are flagged 'unmapped' (these are dropped from the OWL) — typically
    neuprint super-ROIs (SNP(R), INP, NotPrimary, ...), but a large unmapped ROI
    would indicate a missing region mapping."""
    per_roi = counts.groupby("roi")["bodyId"].nunique()
    rows = []
    for roi, neurons in per_roi.items():
        iri, kind = roi2region.get(roi, ("", "unmapped"))
        rows.append((roi, kind, iri, int(neurons)))
    # unmapped ROIs first (they need attention), then alphabetical within group
    rows.sort(key=lambda r: (r[1] != "unmapped", r[0]))
    with open(path, "w") as fh:
        fh.write("roi\tkind\ttarget\tneurons\n")
        for roi, kind, iri, neurons in rows:
            fh.write("%s\t%s\t%s\t%d\n" % (roi, kind, iri, neurons))
    unmapped = [r for r in rows if r[1] == "unmapped"]
    print("ROI mapping: %d ROIs, %d unmapped -> %s"
          % (len(rows), len(unmapped), path))
    if unmapped:
        print("  unmapped:", ", ".join(r[0] for r in unmapped))
    return rows


def main():
    dataset_key, token_file, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    roi_report = sys.argv[4] if len(sys.argv) > 4 else None
    cfg = DATASETS[dataset_key]
    token = open(os.path.expanduser(token_file)).read().strip()

    neuprint_dataset = get_neuprint_dataset(cfg["site"])
    print("neuprint dataset:", neuprint_dataset)

    bodyid2neuron = get_bodyid_to_neuron(cfg["site"])
    print("VFB neurons:", len(bodyid2neuron))

    roi_tsv = cfg.get("roi_mapping_tsv")
    if roi_tsv:
        roi2region = load_roi_mapping_tsv(roi_tsv)
    elif cfg.get("roi_dataset"):
        roi2region = get_roi_to_region(cfg["roi_dataset"])
    else:
        raise RuntimeError(f"No region mapping configured for {dataset_key}")
    print("regions mapped:", len(roi2region))

    counts = fetch_roi_counts(neuprint_dataset, token, bodyid2neuron.keys())
    print("neuprint roi-count rows:", len(counts))

    if roi_report:
        write_roi_report(roi_report, counts, roi2region)

    axioms, used_props = build_axioms(counts, bodyid2neuron, roi2region)
    print("axioms:", len(axioms))

    write_ontology(out_path, dataset_key, neuprint_dataset, axioms, used_props)
    print("wrote", out_path)


if __name__ == "__main__":
    main()

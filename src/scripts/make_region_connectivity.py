"""Generate a region-connectivity OWL file for one connectomics dataset.

For each neuron in a dataset that has a VFB individual, this emits the neuron's
per-neuropil-region synaptic tallies as axiom-annotated relations, matching the
representation already in VFB's production graph (pdb):

  * has presynaptic terminals in  (RO:0002113)  -> annotations: Tbars, downstream
  * has postsynaptic terminal in  (RO:0002110)  -> annotation:  upstream

Counts come from the dataset's source (cfg["source"]): neuprint
(neuprint-python, for the dataset version matching the VFB Site node), the
FlyWire 783 per-neuron neuropil counts (Zenodo), or BANC's per-neuropil
connection counts parquet. flywire and banc have no T-bars, so their
presynaptic edge carries `downstream` (output synapses) instead. The count annotations use the n2o `custom:` namespace
so that `owl2neo4jcsv` turns each axiom into a Neo4j edge carrying the counts as
edge properties (side-load into pdb via vfb-pipeline-dumps).

Each edge is an ObjectPropertyAssertion to a VFB region *individual*; the
individual is side-specific, so laterality is carried by the target.

A second file, synapse_totals_<dataset>.owl, carries each neuron's whole-neuron
synapse totals (total_upstream / total_downstream / total_Tbars) as annotation
assertions on the neuron, which owl2neo4jcsv turns into node properties. These
are the denominators for input-proportion connectivity weights. They cannot be
summed from the region edges (unmapped ROIs are dropped and neuprint ROIs
overlap), and they cannot ride in the region file, which is side-loaded
edge-only.

Datasets marked `totals_only` in dataset_config.yaml (CATMAID) produce only the
totals file, from each skeleton's connectors on the VFB-hosted CATMAID server.
BANC's counts (totals and per region) come from its per-neuropil connection
counts parquet.

Usage:
    python3 make_region_connectivity.py <dataset_key> <neuprint_token_file> <out_dir> <report_dir>
"""
import datetime
import os
import re
import shutil
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import parse_qs, unquote, urlsplit

import fsspec
import pandas as pd
import pyarrow.parquet as pq
import requests
import yaml
import neuprint
from neuprint import NeuronCriteria as NC, fetch_neurons

from cypher_query import query_neo4j, auth, _make_session

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

# source count column -> custom annotation local name, per relation, keyed by the
# dataset's count schema (cfg["counts"], default "neuprint"). Each relation has a
# PRIMARY count that gates edge emission — the edge means the neuron has those
# *terminals* in the region, so an edge is only emitted where the primary count is
# present — plus optional secondary counts.
# (relation, primary source col, primary annotation local, ((secondary col, local), ...))
EDGE_SPECS = {
    # neuprint roi_counts_df: presynaptic edge gated on T-bars (pre); downstream
    # is a secondary count. Matches the existing pdb edge properties.
    "neuprint": (
        (RO_POST, "upstream", "upstream", ()),
        (RO_PRE, "pre", "Tbars", (("downstream", "downstream"),)),
    ),
    # flywire has a single synapse count per direction (no T-bars): input synapses
    # -> upstream (postsynaptic edge), output synapses -> downstream (presynaptic
    # edge, gated on the output count).
    "flywire": (
        (RO_POST, "upstream", "upstream", ()),
        (RO_PRE, "downstream", "downstream", ()),
    ),
    # banc: likewise one synapse count per direction, no T-bars.
    "banc": (
        (RO_POST, "upstream", "upstream", ()),
        (RO_PRE, "downstream", "downstream", ()),
    ),
}
# deterministic ordering of count annotations within an axiom
_LOCAL_ORDER = {"upstream": 0, "Tbars": 1, "downstream": 2}

# Each edge carries only the merged total per count. Region individuals are
# side-specific, so laterality is recoverable downstream from the region's side
# (via part_of in the kb) joined with the neuron's soma side (has_soma_location).

# Datatype of the count annotations. pdb stores these as float
# (neo4j "LIST OF FLOAT", e.g. [545.0]); we match that. Reports cast the counts
# back to integer for readability.
COUNT_DATATYPE = "http://www.w3.org/2001/XMLSchema#float"

# --- Whole-neuron synapse totals ---------------------------------------------
# source column -> node property local name (custom: namespace), per source.
#  * neuprint: the :Neuron node's own totals, which count every synapse on the
#    body whatever the partner (fragments included) and whatever the ROI.
#  * flywire: summed over neuropils from the FlyWire 783 per-neuron neuropil
#    synapse counts (Zenodo 10676866), which count every synapse on the neuron
#    whatever the partner, unproofread fragments included (cleft score >= 50).
TOTAL_SPECS = {
    "neuprint": (("upstream", "total_upstream"),
                 ("downstream", "total_downstream"),
                 ("pre", "total_Tbars")),
    "flywire": (("upstream", "total_upstream"),
                ("downstream", "total_downstream")),
    "banc": (("upstream", "total_upstream"),
             ("downstream", "total_downstream")),
    # catmaid: input sites, and the connectors the skeleton is presynaptic at
    # (CATMAID's T-bar equivalent); synapses made would need every connector's
    # partners fetched, so there is no total_downstream.
    "catmaid": (("upstream", "total_upstream"),
                ("pre", "total_Tbars")),
}
TOTALS_DESCRIPTION = {
    "neuprint": "whole-neuron totals from the neuprint :Neuron node (all synapses, "
                "any partner, any ROI)",
    "flywire": "whole-neuron totals summed from the FlyWire 783 per-neuron neuropil "
               "synapse counts (Zenodo 10676866; all synapses, any partner including "
               "unproofread fragments, cleft score >= 50)",
    "banc": "whole-neuron totals summed from the BANC v1 synapse set per-neuropil "
            "connection counts (synapse size > 3, no connection-count threshold; all "
            "synapses, any partner including unproofread fragments and unsegmented "
            "partners)",
    "catmaid": "whole-neuron totals from the skeleton's connectors on the CATMAID "
               "server (input = postsynaptic links, Tbars = presynaptic connectors; "
               "any partner, traced or not), as of the build date",
}
# Node properties take integers (cf. the scRNAseq cell_count node property);
# only the edge annotations above are float, to match what pdb already holds.
TOTAL_DATATYPE = "http://www.w3.org/2001/XMLSchema#integer"
# concurrent requests to a CATMAID server (one request per skeleton)
CATMAID_WORKERS = 8

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


def _site_link_texts(site, url):
    """The Site node's link_base / homepage strings (source-data links)."""
    df = query_neo4j(
        'MATCH (s:Site {short_form:"%s"}) RETURN s.link_base AS lb, s.homepage AS hp'
        % site,
        url=url,
    )
    if df is False or df.empty:
        raise RuntimeError(f"Site node {site} not found")
    texts = []
    for col in ("lb", "hp"):
        val = df.iloc[0][col]
        text = val[0] if isinstance(val, list) and val else val
        if isinstance(text, str):
            texts.append(unquote(text))  # link_base may be URL-encoded (%3A -> :)
    return texts


def get_neuprint_dataset(site, url=KB):
    """Derive the neuprint dataset string (e.g. 'hemibrain:v1.2.1') from the
    Site node, so the counts come from the version VFB is pinned to."""
    for text in _site_link_texts(site, url):
        m = re.search(r"dataset=([^&\s]+)", text)
        if m:
            return m.group(1)
    raise RuntimeError(f"Could not derive neuprint dataset from Site {site}")


def get_flywire_provenance(site, url=KB, name="flywire"):
    """Provenance string '<name>:<data_version>' from a codex-linked Site node
    (link_base carries data_version=<n>, e.g. flywire 783, banc 888). Used as
    dcterms:source and for the versionInfo, mirroring the neuprint
    '<name>:<version>' shape."""
    for text in _site_link_texts(site, url):
        m = re.search(r"data_version=(\d+)", text)
        if m:
            return "%s:%s" % (name, m.group(1))
    raise RuntimeError(f"Could not derive {name} version from Site {site}")


def get_catmaid_server(site, url=KB):
    """(base URL, project id) of the CATMAID server a Site links to, e.g.
    ('https://l1em.catmaid.virtualflybrain.org', 1)."""
    for text in _site_link_texts(site, url):
        parts = urlsplit(text)
        pid = parse_qs(parts.query).get("pid")
        if "catmaid" in parts.netloc and pid:
            return "%s://%s" % (parts.scheme, parts.netloc), int(pid[0])
    raise RuntimeError(f"Could not derive CATMAID server from Site {site}")


def _local_copy(source, workdir):
    """A local path for source: itself if it is one, else a download into
    workdir. Downloads go through a session that retries 5xx responses, since
    the Zenodo files intermittently 504; pandas' own URL reading does not
    retry."""
    if not re.match(r"https?://", source):
        return source
    path = os.path.join(workdir, os.path.basename(urlsplit(source).path))
    # Zenodo times out about half of requests at its 30 s gateway, so allow
    # more retries than the Neo4j queries need (backoff 2s, 4s, ... capped).
    with _make_session(retries=10, backoff_factor=2) as session, \
            session.get(source, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(path, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
    return path


# the per-neuron neuropil count files label synapses outside every neuropil
# "None"; codex calls them UNASGD, which is kept so ROI reports stay comparable
FLYWIRE_UNASSIGNED = {"None": "UNASGD"}


def load_flywire_counts(post_source, pre_source, body_ids):
    """(totals, counts) from the FlyWire per-neuron neuropil synapse counts
    (feather; one row per neuron x neuropil, columns post_pt_root_id |
    pre_pt_root_id, neuropil, count).

    counts: [bodyId, roi, upstream, downstream], upstream from the post file
    (synapses received in the neuropil), downstream from the pre file
    (synapses made). totals: the same summed over neuropils; every synapse is
    assigned to exactly one neuropil (or UNASGD), so that is the neuron total.
    These count every synapse whatever the partner, unproofread fragments
    included, with no connection-count threshold -- unlike codex's
    neuropil_synapse_table, which only counts connections of >= 5 synapses."""
    ids = {int(b) for b in body_ids}
    workdir = tempfile.mkdtemp(prefix="flywire_counts_")

    def per_neuropil(source, id_col):
        df = pd.read_feather(_local_copy(source, workdir), columns=[id_col, "neuropil", "count"])
        df = df[df[id_col].isin(ids)]
        df = df.assign(neuropil=df["neuropil"].replace(FLYWIRE_UNASSIGNED))
        out = df.groupby([id_col, "neuropil"])["count"].sum()
        out.index.names = ["bodyId", "roi"]
        return out

    try:
        counts = pd.DataFrame({"upstream": per_neuropil(post_source, "post_pt_root_id"),
                               "downstream": per_neuropil(pre_source, "pre_pt_root_id")}
                              ).reset_index()
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    totals = counts.groupby("bodyId", as_index=False)[["upstream", "downstream"]].sum(min_count=1)
    return totals, counts


def load_banc_counts(source, body_ids):
    """(totals, counts) from BANC's per-neuropil connection counts parquet
    (pre_root_id, post_root_id, neuropil, num_synapses; ~177M rows, so read one
    row group at a time).

    totals: whole-neuron [bodyId, upstream, downstream]. counts: per neuron and
    BANC neuropil label [bodyId, neuropil, upstream, downstream]; see
    banc_roi_name for turning the labels into VFB region names.

    This is the v1 synapse set (size > 3) with no connection-count threshold,
    i.e. the synapses VFB's BANC connectivity is built from before its
    count >= 3 cut. Rows whose partner is a fragment or unsegmented (root 0)
    still count: they are synapses on the neuron."""
    ids = pd.Series(sorted(int(b) for b in body_ids), dtype="int64")
    ups, downs = [], []
    with fsspec.open(source, "rb") as fh:
        pf = pq.ParquetFile(fh)
        for i in range(pf.num_row_groups):
            rg = pf.read_row_group(
                i, columns=["pre_root_id", "post_root_id", "neuropil", "num_synapses"]).to_pandas()
            post = rg[rg["post_root_id"].isin(ids)]
            ups.append(post.groupby(["post_root_id", "neuropil"])["num_synapses"].sum())
            pre = rg[rg["pre_root_id"].isin(ids)]
            downs.append(pre.groupby(["pre_root_id", "neuropil"])["num_synapses"].sum())

    def merged(parts):
        if not parts:
            return pd.Series(dtype="int64")
        out = pd.concat(parts).groupby(level=[0, 1]).sum()
        out.index.names = ["bodyId", "neuropil"]
        return out

    counts = pd.DataFrame({"upstream": merged(ups), "downstream": merged(downs)}).reset_index()
    totals = counts.groupby("bodyId", as_index=False)[["upstream", "downstream"]].sum(min_count=1)
    return totals, counts


# BANC v1 neuropil names -> VFB (Bates2026roi) region names. BANC built its v1
# synapse neuropils from the Ito (brain) and Court (VNC) meshes of its
# region_outlines layer, then renamed them into "{region}_{neuropil}_{side}"
# labels (synapses/v1.1/finalize_neuropil_assignments.py in the BANC bucket);
# this inverts those renamings. Court's mVAC was split into T1/T2/T3 by
# position, so all three go back to mVAC; IntTct and LTct exist in both the
# Court and MANC atlases in VFB, and these are the Court ones.
BANC_NEUROPIL_RENAMES = {
    "T1-ProNm": "ProNM-T1", "T2-MesoNm": "MesoNM-T2", "T3-MetaNm": "MetaNM-T3",
    "ANm": "ABDNM", "AMNp": "AMNP",
    "MB-CA": "MB_CA", "MB-ML": "MB_ML", "MB-VL": "MB_VL", "MB-PED": "MB_PED",
    "T1-mVAC": "mVAC", "T2-mVAC": "mVAC", "T3-mVAC": "mVAC",
    "IntTct": "IntTct (COURT)", "LTct": "LTct (COURT)",
}


def banc_roi_name(label, roi2region):
    """VFB region name for a BANC v1 neuropil label like "CB_SLP_R" or
    "VNC_T1-ProNm_L": the sided region (SLP_R) where VFB has one, else the
    unsided one (ProNM-T1), so both sides' synapses land on it. None for the
    unassigned labels and for anything VFB has no region for."""
    _area, rest = label.split("_", 1)
    neuropil, side = rest.rsplit("_", 1)
    name = BANC_NEUROPIL_RENAMES.get(neuropil, neuropil)
    for candidate in ("%s_%s" % (name, side), name):
        if candidate in roi2region:
            return candidate
    return None


def banc_region_counts(counts, roi2region):
    """Per-neuron BANC counts keyed by VFB region name, summing labels that
    resolve to the same region (both sides of an unsided region, mVAC's three
    parts). Labels with no region keep their BANC name so the ROI report lists
    them as unmapped. Also returns each name's brain/vnc group."""
    names = {lab: banc_roi_name(lab, roi2region) or lab for lab in counts["neuropil"].unique()}
    group = {names[lab]: ("vnc" if lab.startswith("VNC_") else "brain") for lab in names}
    out = (counts.assign(roi=counts["neuropil"].map(names))
           .groupby(["bodyId", "roi"], as_index=False)[["upstream", "downstream"]]
           .sum(min_count=1))
    return out, group


def _catmaid_connector_counts(session, base, pid, skid):
    """(upstream, pre) for one skeleton: its postsynaptic links and the
    connectors it is presynaptic at. compact-detail connector rows are
    [treenode, connector, relation, x, y, z] with relation 0 = presynaptic_to,
    1 = postsynaptic_to."""
    r = session.get("%s/%d/skeletons/%d/compact-detail" % (base, pid, skid),
                    params={"with_connectors": "true"}, timeout=120)
    r.raise_for_status()
    rels = [c[2] for c in r.json()[1]]
    return rels.count(1), rels.count(0)


def load_catmaid_totals(base, pid, body_ids, workers=CATMAID_WORKERS):
    """Whole-neuron [bodyId, upstream, pre] totals from a CATMAID server, one
    request per skeleton. A skeleton the server cannot serve (e.g. since
    deleted) is skipped and reported rather than failing the build."""
    session = requests.Session()
    skids = [int(b) for b in body_ids]

    def one(skid):
        try:
            return (skid, *_catmaid_connector_counts(session, base, pid, skid))
        except (requests.RequestException, ValueError, IndexError) as e:
            print("  catmaid skeleton %d skipped: %s" % (skid, e))
            return None

    with ThreadPoolExecutor(max_workers=workers) as ex:
        rows = [r for r in ex.map(one, skids) if r is not None]
    print("catmaid skeletons counted: %d of %d" % (len(rows), len(skids)))
    return pd.DataFrame(rows, columns=["bodyId", "upstream", "pre"])


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
    """Source ROI name -> VFB region-individual IRI.

    A region answers to its label-derived name (see _roi_from_label) and to its
    synonyms: newer ROI datasets record the exact source ROI name as a synonym
    (`synonyms`; hemibrain's use `has_exact_synonym`), which also covers labels
    that differ from it (BANC "LTct (COURT)" has synonym "LTct"). A name that
    resolves to more than one region is ambiguous and dropped, with a warning.
    """
    df = query_neo4j(
        "MATCH (reg)-[:has_source]->(:DataSet {short_form:\"%s\"}) "
        "RETURN reg.iri AS iri, reg.label AS label, "
        "coalesce(reg.synonyms, reg.has_exact_synonym, []) AS syn" % roi_dataset,
        url=url,
    )
    return _roi_names_to_regions(zip(df["iri"], df["label"], df["syn"]))


def _roi_names_to_regions(regions):
    """{name: iri} from (iri, label, synonyms) triples; see get_roi_to_region."""
    names = {}
    for iri, label, syns in regions:
        candidates = {_roi_from_label(label)}
        candidates.update(s.replace("\\", "") for s in (syns or []))
        for name in candidates:
            names.setdefault(name, set()).add(iri)
    ambiguous = sorted(n for n, iris in names.items() if len(iris) > 1)
    if ambiguous:
        print("WARNING: ROI names matching several regions, dropped:", ", ".join(ambiguous))
    return {n: next(iter(iris)) for n, iris in names.items() if len(iris) == 1}


def _roi_from_label(label):
    """VFB region label -> source ROI name.

    Strips the ' on <template>' suffix and any '<Dataset>:' prefix (newer ROI
    datasets label regions e.g. "MANC:LegNp(T3)(R) on JRC2018UnisexVNC"), and
    removes literal backslashes: some VFB labels store an escaped apostrophe
    (e.g. "a\\'L(R)") while neuprint ROI names use a plain apostrophe, so MB
    alpha'/beta' lobe ROIs would otherwise fail to map and be dropped.
    """
    if isinstance(label, list):
        label = label[0]
    name = label.split(" on ")[0].replace("\\", "")
    return re.sub(r"^[A-Za-z0-9]+:", "", name)


NEUPRINT_TOTAL_COLS = ["bodyId", "pre", "post", "downstream", "upstream"]


def fetch_neuprint_counts(client, body_ids):
    """(totals, roi_counts) for the bodyIds we have VFB IRIs for.

    totals is one row per neuron with the :Neuron node's whole-neuron counts;
    roi_counts is neuprint's roi_counts_df (one row per neuron x ROI). Both come
    from the same fetch_neurons call, so they describe the same snapshot."""
    totals, frames = [], []
    ids = list(body_ids)
    chunk = 5000
    for i in range(0, len(ids), chunk):
        ndf, roidf = fetch_neurons(NC(bodyId=ids[i : i + chunk]), client=client)
        if ndf is not None and not ndf.empty:
            ndf = ndf[[c for c in NEUPRINT_TOTAL_COLS if c in ndf.columns]]
            if "upstream" not in ndf.columns and "post" in ndf.columns:
                # older neuprint datasets carry only post (PSD count)
                ndf = ndf.assign(upstream=ndf["post"])
            totals.append(ndf)
        if roidf is not None and not roidf.empty:
            frames.append(roidf)
    totals = (pd.concat(totals, ignore_index=True) if totals
              else pd.DataFrame(columns=NEUPRINT_TOTAL_COLS))
    counts = (pd.concat(frames, ignore_index=True) if frames
              else pd.DataFrame(columns=["bodyId", "roi", "pre", "post", "downstream", "upstream"]))
    return totals, counts


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
    root = client.meta["roiHierarchy"]
    groups = roi_groups_from_hierarchy(root)
    # Some roiInfo ROIs are absent from the hierarchy tree (e.g. hemibrain
    # FB-column3). They can't be under the VNC branch, so default them to the
    # root's group (brain, unless the whole dataset is VNC-rooted like manc).
    root_group = groups.get(root.get("name"), "brain")
    for roi in client.meta.get("roiInfo", {}):
        groups.setdefault(roi, root_group)
    return groups


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


def build_axioms(counts, bodyid2neuron, roi2region, edge_specs=EDGE_SPECS["neuprint"]):
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
        for rel, pcol, plocal, secondary in edge_specs:
            cs = _edge_counts(row, pcol, plocal, secondary)
            if cs is not None:
                axioms.append(_relation_axiom(neuron, target, rel, _count_annotations(cs)))
                used_props.add(rel)
    return axioms, used_props


def build_total_axioms(totals, bodyid2neuron, total_spec=TOTAL_SPECS["neuprint"]):
    """One AnnotationAssertion per (neuron, total) on the neuron individual, so
    owl2neo4jcsv makes it a node property. Zero/missing totals are omitted, as
    are neurons with no VFB individual."""
    axioms = []
    # Column-wise, not iterrows(): a row holding a NaN is upcast to float64,
    # which corrupts 18-digit flywire/banc root ids so they no longer match.
    ids = totals["bodyId"].astype("int64").tolist()
    cols = [(totals[c].tolist() if c in totals.columns else [None] * len(ids), local)
            for c, local in total_spec]
    for i, body in enumerate(ids):
        neuron = bodyid2neuron.get(body)
        if neuron is None:
            continue
        for values, local in cols:
            v = values[i]
            if pd.notna(v) and v > 0:
                axioms.append('AnnotationAssertion(<%s%s> <%s> "%d"^^<%s>)'
                              % (CUSTOM, local, neuron, int(v), TOTAL_DATATYPE))
    return axioms


REGION_ONTOLOGY = ("region_connectivity",
                   "VFB region connectivity: %s",
                   "Region (neuropil) synaptic connectivity for %s, generated from %s.")
TOTALS_ONTOLOGY = ("synapse_totals",
                   "VFB neuron synapse totals: %s",
                   "Whole-neuron synapse totals for %s, generated from %s.")


def write_ontology(path, dataset_key, source_id, axioms, used_props,
                   kind=REGION_ONTOLOGY, extra_comment=None):
    name, title, comment = kind
    onto_iri = "http://virtualflybrain.org/data/VFB/OWL/%s_%s.owl" % (name, dataset_key)
    comment = comment % (dataset_key, source_id)
    if extra_comment:
        comment += " Counts are " + extra_comment + "."
    # declare exactly the custom annotation properties actually used
    used_anns = set(re.findall(r"<http://n2o\.neo/custom/([A-Za-z_]+)>", "\n".join(axioms)))
    decls = [
        "Declaration(AnnotationProperty(<%s%s>))" % (CUSTOM, n)
        for n in sorted(used_anns)
    ]
    decls += ["Declaration(ObjectProperty(<%s>))" % p for p in sorted(used_props)]
    # source_id is "<name>:<version>", e.g. hemibrain:v1.2.1 or flywire:783
    version = source_id.partition(":")[2] or source_id
    lines = [
        "Prefix(:=<http://purl.obolibrary.org/obo/>)",
        "Prefix(owl:=<http://www.w3.org/2002/07/owl#>)",
        "Prefix(xsd:=<http://www.w3.org/2001/XMLSchema#>)",
        "Prefix(dcterms:=<http://purl.org/dc/terms/>)",
        "Prefix(custom:=<%s>)" % CUSTOM,
        "",
        "Ontology(<%s>" % onto_iri,
        'Annotation(rdfs:comment "%s")' % comment,
        'Annotation(dcterms:title "%s")' % (title % dataset_key),
        'Annotation(dcterms:identifier "%s")' % dataset_key,
        'Annotation(dcterms:source "%s")' % source_id,
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


def report_totals(totals, counts, roi2region):
    """Print how the whole-neuron input totals compare with the summed mapped
    region input counts. Where the mapped regions partition the neuron
    (flywire, banc) the region sum can only be smaller, by the synapses in
    unmapped regions; where they overlap (hemibrain and male-cns map super-ROIs
    as well as their parts) it can be larger. A large unexplained shift between
    builds is the thing to look for."""
    if totals.empty:
        print("WARNING: no whole-neuron totals fetched")
        return
    region_up = (counts[counts["roi"].isin(roi2region)]
                 .groupby("bodyId")["upstream"].sum())
    both = totals.set_index("bodyId")["upstream"].to_frame("total").join(
        region_up.rename("regions"), how="inner")
    print("neurons with totals: %d (with region counts too: %d)" % (len(totals), len(both)))
    if len(both):
        print("  total_upstream sum %d vs summed mapped-region upstream %d; "
              "%d neurons have total < region sum"
              % (both["total"].sum(), both["regions"].sum(),
                 int((both["total"] < both["regions"]).sum())))


def main():
    # args: <dataset> <token_file> <out_dir> <report_dir>
    dataset_key, token_file, out_dir, report_dir = sys.argv[1:5]
    if dataset_key not in DATASETS:
        sys.exit("unknown dataset %r; configured: %s"
                 % (dataset_key, ", ".join(sorted(DATASETS))))
    cfg = DATASETS[dataset_key]
    source = cfg.get("source", "neuprint")

    # Default to the kb: dataset updates land there first (needed to migrate to
    # a new dataset version before it reaches pdb).
    url = ENDPOINTS[cfg.get("endpoint", "kb")]

    bodyid2neuron = get_bodyid_to_neuron(cfg["site"], url)
    print("VFB neurons:", len(bodyid2neuron))
    # totals-only datasets (no region counts loaded for them yet) skip the
    # region file, its ROI mapping and its report
    totals_only = cfg.get("totals_only", False)
    roi2region = {} if totals_only else get_roi_to_region(cfg["roi_dataset"], url)
    if not totals_only:
        print("regions mapped:", len(roi2region))

    counts = None
    if source == "neuprint":
        source_id = get_neuprint_dataset(cfg["site"], url)
        print("neuprint dataset:", source_id)
        token = open(os.path.expanduser(token_file)).read().strip()
        client = neuprint.Client("https://neuprint.janelia.org", dataset=source_id, token=token)
        roi_group = get_roi_groups(client)
        totals, counts = fetch_neuprint_counts(client, bodyid2neuron.keys())
    elif source == "flywire":
        source_id = get_flywire_provenance(cfg["site"], url)
        print("flywire source:", source_id)
        totals, counts = load_flywire_counts(cfg["counts_post_url"], cfg["counts_pre_url"],
                                             bodyid2neuron.keys())
        roi_group = {r: "brain" for r in roi2region}  # FAFB is brain-only (no VNC)
    elif source == "banc":
        source_id = get_flywire_provenance(cfg["site"], url, name="banc")
        print("banc source:", source_id)
        totals, banc_counts = load_banc_counts(cfg["counts_parquet_url"], bodyid2neuron.keys())
        if not totals_only:
            counts, roi_group = banc_region_counts(banc_counts, roi2region)
    elif source == "catmaid" and totals_only:
        base, pid = get_catmaid_server(cfg["site"], url)
        # a live server has no version: the build date stands in for one
        source_id = "%s:%s" % (cfg["site"], datetime.date.today().isoformat())
        print("catmaid source:", base, "pid", pid, "->", source_id)
        totals = load_catmaid_totals(base, pid, bodyid2neuron.keys())
    else:
        sys.exit("unsupported source %r for dataset %s (catmaid needs totals_only)"
                 % (source, dataset_key))

    qc_dir = os.path.join(out_dir, "qc")
    os.makedirs(qc_dir, exist_ok=True)
    if not totals_only:
        print("region-count rows:", len(counts))
        os.makedirs(report_dir, exist_ok=True)
        # roi_group only annotates the report's `group` column (brain/vnc); one file
        # per dataset (a whole-CNS dataset stays well under GitHub's size limit).
        write_roi_report(os.path.join(report_dir, "roi_mapping.tsv"), counts, roi2region, roi_group)

        axioms, used_props = build_axioms(counts, bodyid2neuron, roi2region, EDGE_SPECS[source])
        print("axioms:", len(axioms))

        out_path = os.path.join(out_dir, "region_connectivity_%s.ofn" % dataset_key)
        write_ontology(out_path, dataset_key, source_id, axioms, used_props)
        # small sample for build-time QC (full file is too large for robot verify)
        write_ontology(os.path.join(qc_dir, "region_connectivity_%s.ofn" % dataset_key),
                       dataset_key, source_id, axioms[:SAMPLE_AXIOMS], used_props)
        print("wrote", out_path)

    # Whole-neuron totals: a separate file because it is loaded into the main
    # graph (node properties), whereas the region file is side-loaded edge-only.
    total_axioms = build_total_axioms(totals, bodyid2neuron, TOTAL_SPECS[source])
    if counts is not None:
        report_totals(totals, counts, roi2region)
    else:
        print("neurons with totals: %d" % len(totals))
    totals_path = os.path.join(out_dir, "synapse_totals_%s.ofn" % dataset_key)
    for path, axs in ((totals_path, total_axioms),
                      (os.path.join(qc_dir, "synapse_totals_%s.ofn" % dataset_key),
                       total_axioms[:SAMPLE_AXIOMS])):
        write_ontology(path, dataset_key, source_id, axs, set(),
                       kind=TOTALS_ONTOLOGY, extra_comment=TOTALS_DESCRIPTION[source])
    print("total axioms:", len(total_axioms))
    print("wrote", totals_path)

if __name__ == "__main__":
    main()

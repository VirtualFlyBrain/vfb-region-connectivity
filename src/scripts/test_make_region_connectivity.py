"""Tests for make_region_connectivity.

Unit tests exercise the pure logic; the integration tests run the whole pipeline
on a *saved* neuprint/pdb sample (test_data/) with no live queries.

Run: pip install -r requirements.txt pytest && pytest src/scripts
"""
import json
import os
import shutil
import subprocess
import sys

import pandas as pd
import pytest

import make_region_connectivity as m

CUSTOM = "http://n2o.neo/custom/"
FIX = os.path.join(os.path.dirname(__file__), "test_data")
SPARQL = os.path.join(os.path.dirname(__file__), "..", "sparql")
RC_CHECKS = [
    "post-edge-missing-upstream",
    "pre-edge-missing-count",
    "nonpositive-count",
    "unexpected-connectivity-relation",
    "unexpected-count-annotation",
    "nonvfb-neuron",
    "nonvfb-region",
]


# --- helpers ---------------------------------------------------------------
def test_lit_is_float_typed_integer_string():
    # counts are whole numbers but typed xsd:float to match pdb
    assert m._lit(23.0) == '"23"^^<http://www.w3.org/2001/XMLSchema#float>'
    assert m._lit(1665) == '"1665"^^<http://www.w3.org/2001/XMLSchema#float>'


# --- dataset config --------------------------------------------------------
def test_datasets_loaded_from_yaml():
    cfg = m.load_datasets()
    assert cfg == m.DATASETS                      # module-level load matches
    assert "hemibrain" in cfg and "male-cns" in cfg
    for entry in cfg.values():                    # required keys present
        assert entry["site"]
        assert entry.get("totals_only") or entry["roi_dataset"]   # region datasets need ROIs


# --- ROI grouping (brain/vnc) ----------------------------------------------
def test_roi_groups_root_is_vnc():
    # manc: the hierarchy root itself is 'VNC' -> every ROI is vnc
    h = {"name": "VNC", "children": [{"name": "LTct"}, {"name": "ANm"}]}
    assert m.roi_groups_from_hierarchy(h) == {"VNC": "vnc", "LTct": "vnc", "ANm": "vnc"}


def test_roi_groups_mixed_cns():
    # male-cns/optic-lobe: a 'VNC' branch under a CNS root; rest is brain
    h = {"name": "CNS", "children": [
        {"name": "Brain", "children": [{"name": "ME(R)"}]},
        {"name": "VNC", "children": [{"name": "LegNp(T1)(L)"}]},
    ]}
    g = m.roi_groups_from_hierarchy(h)
    assert g["ME(R)"] == "brain" and g["VNC"] == "vnc" and g["LegNp(T1)(L)"] == "vnc"


def test_roi_groups_brain_only():
    h = {"name": "hemibrain", "children": [{"name": "SIP(R)"}]}
    assert m.roi_groups_from_hierarchy(h) == {"hemibrain": "brain", "SIP(R)": "brain"}


class _StubClient:
    def __init__(self, meta):
        self.meta = meta


def test_get_roi_groups_defaults_rois_absent_from_hierarchy():
    # FB-column3 is in roiInfo but not in the hierarchy tree -> default to root group (brain)
    meta = {"roiHierarchy": {"name": "hemibrain", "children": [{"name": "FB"}]},
            "roiInfo": {"FB": {}, "FB-column3": {}}}
    g = m.get_roi_groups(_StubClient(meta))
    assert g["FB"] == "brain" and g["FB-column3"] == "brain"


def test_get_roi_groups_vnc_rooted_defaults_to_vnc():
    meta = {"roiHierarchy": {"name": "VNC", "children": [{"name": "LTct"}]},
            "roiInfo": {"LTct": {}, "SomeExtra": {}}}
    g = m.get_roi_groups(_StubClient(meta))
    assert g["LTct"] == "vnc" and g["SomeExtra"] == "vnc"


# --- ROI label parsing (regression for the backslash bug) ------------------
def test_roi_from_label_strips_suffix():
    assert m._roi_from_label("SIP(R) on JRC_FlyEM_Hemibrain") == "SIP(R)"


def test_roi_from_label_strips_backslash_escaped_apostrophe():
    # VFB stores a\'L(R); neuprint uses a'L(R)
    assert m._roi_from_label("a\\'L(R) on JRC_FlyEM_Hemibrain") == "a'L(R)"


def test_roi_from_label_handles_list():
    assert m._roi_from_label(["LH(R) on X"]) == "LH(R)"


# --- edge counts (primary gates the edge) ----------------------------------
POST, PRE = m.EDGE_SPECS["neuprint"]
FW_POST, FW_PRE = m.EDGE_SPECS["flywire"]


def test_edge_counts_postsynaptic():
    assert m._edge_counts({"upstream": 5}, *POST[1:]) == {"upstream": 5}
    assert m._edge_counts({"upstream": 0}, *POST[1:]) is None


def test_edge_counts_presynaptic_requires_tbars():
    # T-bars present -> edge; downstream included only when > 0
    assert m._edge_counts({"pre": 3, "downstream": 7}, *PRE[1:]) == {"Tbars": 3, "downstream": 7}
    assert m._edge_counts({"pre": 3, "downstream": 0}, *PRE[1:]) == {"Tbars": 3}
    # no T-bars -> no presynaptic edge even if downstream > 0 (the male-cns case)
    assert m._edge_counts({"pre": 0, "downstream": 4}, *PRE[1:]) is None


def test_edge_counts_flywire_presynaptic_uses_downstream():
    # flywire: no Tbars; the presynaptic edge is gated on the output count (downstream)
    assert m._edge_counts({"upstream": 5, "downstream": 8}, *FW_PRE[1:]) == {"downstream": 8}
    assert m._edge_counts({"upstream": 5, "downstream": 0}, *FW_PRE[1:]) is None
    assert m._edge_counts({"upstream": 5}, *FW_POST[1:]) == {"upstream": 5}


# --- axiom shape -----------------------------------------------------------
def test_relation_axiom_is_object_property_assertion():
    ax = m._relation_axiom("n", "r", m.RO_POST, ['Annotation(x)'])
    assert ax.startswith("ObjectPropertyAssertion(")
    assert "<n>" in ax and "<r>" in ax


# --- build_axioms ----------------------------------------------------------
def _counts():
    return pd.DataFrame([
        {"bodyId": 1, "roi": "SIP(R)", "upstream": 5, "pre": 2, "downstream": 10},
        {"bodyId": 1, "roi": "SCL(R)", "upstream": 3, "pre": 0, "downstream": 0},
        {"bodyId": 1, "roi": "SNP(R)", "upstream": 9, "pre": 9, "downstream": 9},  # unmapped
        {"bodyId": 99, "roi": "SIP(R)", "upstream": 1, "pre": 1, "downstream": 1},  # unknown neuron
    ])


_MAP = {"SIP(R)": "http://vfb/sip", "SCL(R)": "http://vfb/scl"}


def test_build_axioms_counts_and_filtering():
    axioms, props = m.build_axioms(_counts(), {1: "http://vfb/n1"}, _MAP)
    # (1,SIP): post+pre = 2 ; (1,SCL): post only = 1 ; unmapped + unknown neuron dropped
    assert len(axioms) == 3
    assert props == {m.RO_POST, m.RO_PRE}


def test_no_side_buckets():
    axioms, _ = m.build_axioms(_counts(), {1: "http://vfb/n1"}, _MAP)
    assert not any("_left" in a or "_right" in a for a in axioms)


def test_build_axioms_scl_emits_only_postsynaptic():
    axioms, _ = m.build_axioms(_counts(), {1: "http://vfb/n1"}, _MAP)
    scl = [a for a in axioms if "http://vfb/scl" in a]
    assert len(scl) == 1
    assert m.RO_POST in scl[0] and m.RO_PRE not in scl[0]


# --- flywire (CSV source, upstream+downstream, no Tbars) -------------------
def _fw_counts():
    return pd.DataFrame([
        {"bodyId": 1, "roi": "ME_L", "upstream": 12, "downstream": 8},
        {"bodyId": 1, "roi": "EB", "upstream": 0, "downstream": 5},   # pre-only
        {"bodyId": 1, "roi": "LO_R", "upstream": 4, "downstream": 0},  # post-only
    ])


def test_build_axioms_flywire_uses_downstream_no_tbars():
    axioms, props = m.build_axioms(_fw_counts(), {1: "http://vfb/n1"},
                                   {"ME_L": "http://vfb/me", "EB": "http://vfb/eb",
                                    "LO_R": "http://vfb/lo"}, m.EDGE_SPECS["flywire"])
    # ME_L: post+pre=2 ; EB: pre only=1 ; LO_R: post only=1
    assert len(axioms) == 4
    assert props == {m.RO_POST, m.RO_PRE}
    assert not any("Tbars" in a for a in axioms)            # flywire has no T-bars
    pre_axioms = [a for a in axioms if m.RO_PRE in a]
    assert pre_axioms and all("custom/downstream" in a for a in pre_axioms)


def test_load_flywire_totals_sums_over_neuropils(tmp_path):
    post, pre = tmp_path / "post.feather", tmp_path / "pre.feather"
    pd.DataFrame({"post_pt_root_id": [1, 1, 2, 3],
                  "neuropil": ["ME_L", "LO_L", "EB", "EB"],
                  "count": [1, 4, 6, 9]}).to_feather(post)
    pd.DataFrame({"pre_pt_root_id": [1, 2],          # 2 has outputs, 3 none listed
                  "neuropil": ["ME_L", "EB"],
                  "count": [3, 2]}).to_feather(pre)
    df = m.load_flywire_totals(str(post), str(pre), body_ids=[1, 2]).set_index("bodyId")
    assert df.loc[1, "upstream"] == 5 and df.loc[1, "downstream"] == 3
    assert df.loc[2, "upstream"] == 6 and df.loc[2, "downstream"] == 2
    assert 3 not in df.index                          # not requested


def test_load_banc_totals_sums_over_row_groups(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    pqf = tmp_path / "banc.parquet"
    tbl = pa.table({"pre_root_id": [0, 1, 1, 2, 3],      # 0 = unsegmented partner
                    "post_root_id": [1, 2, 2, 1, 1],
                    "neuropil": ["CB_AL_L", "CB_AL_L", "CB_SLP_R", "VNC_x", "CB_AL_L"],
                    "num_synapses": [1, 2, 3, 4, 5]})
    pq.write_table(tbl, pqf, row_group_size=2)            # 3 row groups to merge
    df = m.load_banc_totals(str(pqf), body_ids=[1, 2]).set_index("bodyId")
    assert df.loc[1, "upstream"] == 1 + 4 + 5 and df.loc[1, "downstream"] == 5
    assert df.loc[2, "upstream"] == 5 and df.loc[2, "downstream"] == 4
    assert 3 not in df.index


def test_load_catmaid_totals_counts_links_and_skips_failures(monkeypatch):
    import requests
    fake = {11: (7, 2), 12: (0, 3)}

    def counts(session, base, pid, skid):
        if skid not in fake:
            raise requests.HTTPError("404")
        return fake[skid]
    monkeypatch.setattr(m, "_catmaid_connector_counts", counts)
    df = m.load_catmaid_totals("https://x", 1, [11, 12, 13], workers=2).set_index("bodyId")
    assert df.loc[11, "upstream"] == 7 and df.loc[11, "pre"] == 2
    assert df.loc[12, "upstream"] == 0 and 13 not in df.index   # failed skeleton skipped


def test_catmaid_connector_counts_reads_relations():
    class _Resp:
        def raise_for_status(self): pass
        def json(self):   # [nodes, connectors]; relation 0 = pre, 1 = post
            return [[], [[1, 100, 1, 0, 0, 0], [2, 101, 1, 0, 0, 0], [3, 102, 0, 0, 0, 0]]]

    class _Session:
        def get(self, url, params=None, timeout=None):
            assert url == "https://x/1/skeletons/5/compact-detail"
            return _Resp()
    assert m._catmaid_connector_counts(_Session(), "https://x", 1, 5) == (2, 1)


def test_get_catmaid_server_from_site_link(monkeypatch):
    monkeypatch.setattr(m, "_site_link_texts", lambda site, url: [
        "https://l1em.catmaid.virtualflybrain.org/?pid=1&tool=tracingtool&active_skeleton_id="])
    assert m.get_catmaid_server("catmaid_l1em") == ("https://l1em.catmaid.virtualflybrain.org", 1)


def test_load_flywire_counts_reshapes_wide_table(tmp_path):
    csv = tmp_path / "syn.csv"
    csv.write_text(
        "root_id,input synapses,output synapses,"
        "input synapses in ME_L,output synapses in ME_L,"
        "input synapses in EB,output synapses in EB,"
        "input partners in ME_L,output partners in ME_L\n"
        "1,16,5,12,0,4,5,3,0\n"        # ME_L post-only, EB pre-only
        "2,0,0,0,0,0,0,0,0\n"          # neuron with no synapses -> no rows
        "3,1,1,0,0,0,0,0,0\n")         # root_id 3 not requested -> excluded
    df = m.load_flywire_counts(str(csv), body_ids=[1, 2])
    rows = {(int(r.bodyId), r.roi): (r.upstream, r.downstream) for r in df.itertuples()}
    assert rows == {(1, "ME_L"): (12, 0), (1, "EB"): (4, 5)}   # partners ignored
    assert 3 not in {b for b, _ in rows}


# --- whole-neuron synapse totals --------------------------------------------
def test_build_total_axioms_neuprint():
    totals = pd.DataFrame([
        {"bodyId": 1, "pre": 3, "post": 40, "upstream": 40, "downstream": 12},
        {"bodyId": 2, "pre": 0, "post": 5, "upstream": 5, "downstream": 0},
        {"bodyId": 99, "pre": 1, "post": 1, "upstream": 1, "downstream": 1},  # unknown neuron
    ])
    axioms = m.build_total_axioms(totals, {1: "http://vfb/n1", 2: "http://vfb/n2"})
    assert len(axioms) == 4                       # n1: 3 totals ; n2: upstream only
    assert ('AnnotationAssertion(<%stotal_upstream> <http://vfb/n1> '
            '"40"^^<http://www.w3.org/2001/XMLSchema#integer>)' % CUSTOM) in axioms
    assert not any("<http://vfb/n2>" in a and "total_downstream" in a for a in axioms)


def test_build_total_axioms_flywire_has_no_tbars():
    totals = pd.DataFrame([{"bodyId": 1, "upstream": 9, "downstream": 5}])
    axioms = m.build_total_axioms(totals, {1: "http://vfb/n1"}, m.TOTAL_SPECS["flywire"])
    assert len(axioms) == 2 and not any("Tbars" in a for a in axioms)


def test_build_total_axioms_keeps_large_root_ids_with_nan():
    # a NaN in the row must not turn an 18-digit root id into a lossy float
    rid = 720575940618201620
    totals = pd.DataFrame({"bodyId": [rid], "upstream": [150.0], "downstream": [float("nan")]})
    axioms = m.build_total_axioms(totals, {rid: "http://vfb/n1"}, m.TOTAL_SPECS["flywire"])
    assert len(axioms) == 1 and "total_upstream" in axioms[0]


def test_fetch_neuprint_counts_falls_back_to_post(monkeypatch):
    # a dataset whose :Neuron has no upstream property: post stands in for it
    ndf = pd.DataFrame([{"bodyId": 1, "pre": 2, "post": 30, "downstream": 8, "instance": "x"}])
    roidf = pd.DataFrame([{"bodyId": 1, "roi": "SIP(R)", "pre": 2, "post": 30,
                           "downstream": 8, "upstream": 30}])
    monkeypatch.setattr(m, "fetch_neurons", lambda crit, client=None: (ndf, roidf))
    totals, counts = m.fetch_neuprint_counts(object(), [1])
    assert totals.loc[0, "upstream"] == 30 and "instance" not in totals.columns
    assert len(counts) == 1


# --- reports & ontology ----------------------------------------------------
def test_write_roi_report_flags_unmapped(tmp_path):
    out = tmp_path / "roi_mapping.tsv"
    m.write_roi_report(str(out), _counts(), _MAP, {"SIP(R)": "brain", "SNP(R)": "brain"})
    lines = out.read_text().splitlines()
    assert lines[0] == "roi\tgroup\ttarget\tneurons"
    body = {ln.split("\t")[0]: ln.split("\t") for ln in lines[1:]}
    assert body["SNP(R)"][2] == "" and body["SNP(R)"][1] == "unmapped"   # unmapped: empty target, group flags it
    assert body["SIP(R)"][2] == "http://vfb/sip" and body["SIP(R)"][1] == "brain"
    # unmapped ROIs float to the top
    assert lines[1].split("\t")[2] == ""


def test_write_ontology_header_has_dataset_and_version(tmp_path):
    out = tmp_path / "x.ofn"
    m.write_ontology(str(out), "hemibrain", "hemibrain:v1.2.1", ["ClassAssertion(a b)"], {m.RO_POST})
    text = out.read_text()
    assert 'dcterms:identifier "hemibrain"' in text
    assert 'dcterms:source "hemibrain:v1.2.1"' in text
    assert 'owl:versionInfo "v1.2.1"' in text
    assert "ClassAssertion(a b)" in text


def test_write_totals_ontology_header(tmp_path):
    out = tmp_path / "t.ofn"
    m.write_ontology(str(out), "flywire", "flywire:783", [], set(),
                     kind=m.TOTALS_ONTOLOGY, extra_comment=m.TOTALS_DESCRIPTION["flywire"])
    text = out.read_text()
    assert "<http://virtualflybrain.org/data/VFB/OWL/synapse_totals_flywire.owl>" in text
    assert "Zenodo 10676866" in text               # what the counts mean is recorded


# --- integration: whole pipeline on saved samples (no live queries) ---------
def _run_main(dataset, tmp_path, monkeypatch, counts, b2n, np_dataset, roi_group, r2r):
    """Run main() fully offline: every network call is stubbed with saved data.
    Returns (report_dir, [produced OFN paths])."""
    monkeypatch.setattr(m, "get_neuprint_dataset", lambda site, url=None: np_dataset)
    monkeypatch.setattr(m, "get_bodyid_to_neuron", lambda site, url=None: b2n)
    monkeypatch.setattr(m, "get_roi_to_region", lambda ds, url=None: r2r)
    totals = (counts.groupby("bodyId", as_index=False)[["pre", "post", "upstream", "downstream"]]
              .sum())
    monkeypatch.setattr(m, "fetch_neuprint_counts", lambda client, ids: (totals, counts))
    monkeypatch.setattr(m, "get_roi_groups", lambda client: roi_group)
    monkeypatch.setattr(m.neuprint, "Client", lambda *a, **k: object())  # no connection
    token = tmp_path / "token.txt"
    token.write_text("dummy")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    report_dir = tmp_path / "report"
    monkeypatch.setattr(sys, "argv", ["prog", dataset, str(token), str(out_dir), str(report_dir)])
    m.main()
    return report_dir, sorted(out_dir.glob("region_connectivity_*.ofn"))


def _load_region_map(name):  # fixtures store [iri, kind]; the pipeline uses plain IRIs
    return {k: (v[0] if isinstance(v, list) else v)
            for k, v in json.load(open(os.path.join(FIX, name))).items()}


def _hemibrain(monkeypatch, tmp_path):
    counts = pd.read_csv(os.path.join(FIX, "neuprint_roi_counts_sample.tsv"), sep="\t")
    b2n = {int(k): v for k, v in json.load(open(os.path.join(FIX, "bodyid_to_neuron.json"))).items()}
    r2r = _load_region_map("roi_to_region.json")
    groups = {r: "brain" for r in counts["roi"].unique()}   # hemibrain is brain-only
    report_dir, ofns = _run_main("hemibrain", tmp_path, monkeypatch, counts, b2n,
                                 "hemibrain:v1.2.1", groups, r2r)
    return counts, r2r, report_dir, ofns


def _male_cns(monkeypatch, tmp_path):
    counts = pd.read_csv(os.path.join(FIX, "male-cns_roi_counts_sample.tsv"), sep="\t")
    b2n = {int(k): v for k, v in json.load(open(os.path.join(FIX, "male-cns_bodyid_to_neuron.json"))).items()}
    groups = json.load(open(os.path.join(FIX, "male-cns_roi_groups.json")))
    r2r = _load_region_map("male-cns_roi_to_region.json")   # Berg2025aroi individuals (per-side)
    return _run_main("male-cns", tmp_path, monkeypatch, counts, b2n, "male-cns:v1.0", groups, r2r)


def test_male_cns_individual_path(tmp_path, monkeypatch):
    _report_dir, ofns = _male_cns(monkeypatch, tmp_path)
    # one file per dataset (brain + vnc regions together)
    assert [o.name for o in ofns] == ["region_connectivity_male-cns.ofn"]
    text = ofns[0].read_text()
    assert "ObjectPropertyAssertion(" in text          # individual (per-side) edges
    assert "_left" not in text and "_right" not in text  # no side-specific annotations


def _flywire(monkeypatch, tmp_path):
    """flywire main() offline: kb lookups + the CSV load are stubbed."""
    b2n = {1: "http://virtualflybrain.org/reports/VFB_fw00000001"}
    r2r = {"ME_L": "http://virtualflybrain.org/reports/VFB_fw000000me",
           "EB": "http://virtualflybrain.org/reports/VFB_fw000000eb",
           "LO_R": "http://virtualflybrain.org/reports/VFB_fw000000lo"}
    monkeypatch.setattr(m, "get_flywire_provenance", lambda site, url=None: "flywire:783")
    monkeypatch.setattr(m, "get_bodyid_to_neuron", lambda site, url=None: b2n)
    monkeypatch.setattr(m, "get_roi_to_region", lambda ds, url=None: r2r)
    monkeypatch.setattr(m, "load_flywire_counts", lambda src, ids: _fw_counts())
    monkeypatch.setattr(m, "load_flywire_totals", lambda post, pre, ids: pd.DataFrame(
        [{"bodyId": 1, "upstream": 30, "downstream": 20}]))
    token = tmp_path / "token.txt"; token.write_text("dummy")
    out_dir = tmp_path / "out"; out_dir.mkdir()
    report_dir = tmp_path / "report"
    monkeypatch.setattr(sys, "argv", ["prog", "flywire", str(token), str(out_dir), str(report_dir)])
    m.main()
    return report_dir, sorted(out_dir.glob("region_connectivity_*.ofn"))


def test_totals_only_dataset_writes_no_region_file(tmp_path, monkeypatch):
    monkeypatch.setitem(m.DATASETS, "catmaid-test",
                        {"site": "catmaid_l1em", "source": "catmaid", "totals_only": True})
    monkeypatch.setattr(m, "get_bodyid_to_neuron",
                        lambda site, url=None: {11: "http://virtualflybrain.org/reports/VFB_00000011"})
    monkeypatch.setattr(m, "get_roi_to_region", lambda *a, **k: pytest.fail("no ROIs needed"))
    monkeypatch.setattr(m, "get_catmaid_server", lambda site, url=None: ("https://x", 1))
    monkeypatch.setattr(m, "load_catmaid_totals", lambda base, pid, ids: pd.DataFrame(
        [{"bodyId": 11, "upstream": 7, "pre": 2}]))
    out_dir = tmp_path / "out"; out_dir.mkdir()
    monkeypatch.setattr(sys, "argv", ["prog", "catmaid-test", "unused", str(out_dir),
                                      str(tmp_path / "report")])
    m.main()
    assert not list(out_dir.glob("region_connectivity_*"))
    text = (out_dir / "synapse_totals_catmaid-test.ofn").read_text()
    assert "total_upstream>" in text and "total_Tbars>" in text and "total_downstream>" not in text
    assert 'dcterms:source "catmaid_l1em:' in text                # build date stands in for a version


def test_flywire_path(tmp_path, monkeypatch):
    report_dir, ofns = _flywire(monkeypatch, tmp_path)
    assert [o.name for o in ofns] == ["region_connectivity_flywire.ofn"]
    text = ofns[0].read_text()
    assert 'dcterms:source "flywire:783"' in text and 'owl:versionInfo "783"' in text
    assert "ObjectPropertyAssertion(" in text
    assert "custom/downstream" in text and "custom/Tbars" not in text  # no T-bars
    # FAFB is brain-only: mapped ROIs grouped 'brain'
    rows = {ln.split("\t")[0]: ln.split("\t")
            for ln in (report_dir / "roi_mapping.tsv").read_text().splitlines()[1:]}
    assert rows["ME_L"][1] == "brain"
    totals = (ofns[0].parent / "synapse_totals_flywire.ofn").read_text()
    assert totals.count("AnnotationAssertion(") == 2 and "total_Tbars" not in totals


def test_individual_pipeline_from_sample(tmp_path, monkeypatch):
    counts, r2r, report_dir, ofns = _hemibrain(monkeypatch, tmp_path)
    # brain-only -> a single, unsuffixed file
    assert [o.name for o in ofns] == ["region_connectivity_hemibrain.ofn"]
    text = ofns[0].read_text()
    assert 'dcterms:identifier "hemibrain"' in text and 'owl:versionInfo "v1.2.1"' in text
    mapped = counts[counts["roi"].isin(r2r)]
    assert text.count("<%s>" % m.RO_POST) == int((mapped["upstream"] > 0).sum()) + 1
    assert text.count("<%s>" % m.RO_PRE) == int((mapped["pre"] > 0).sum()) + 1
    rows = {ln.split("\t")[0]: ln.split("\t")
            for ln in (report_dir / "roi_mapping.tsv").read_text().splitlines()[1:]}
    for sup in ("INP", "NotPrimary", "OL(R)"):
        assert rows[sup][2] == ""   # unmapped super-ROIs: empty target
    # one totals file alongside, with an upstream total for every known neuron
    totals = (ofns[0].parent / "synapse_totals_hemibrain.ofn").read_text()
    known = set(counts["bodyId"]) & set(int(b) for b in json.load(open(os.path.join(FIX, "bodyid_to_neuron.json"))))
    assert totals.count("AnnotationAssertion(<%stotal_upstream>" % CUSTOM) == len(known)


def _shape_check(ofn, tmp_path):
    """Convert an OFN and run the shape checks; return robot verify's exit code."""
    owl = tmp_path / (ofn.stem + ".owl")
    subprocess.run(["robot", "convert", "-i", str(ofn), "-f", "owl", "-o", str(owl)], check=True)
    queries = [os.path.join(SPARQL, "%s-violation.sparql" % c) for c in RC_CHECKS]
    return subprocess.run(
        ["robot", "verify", "-i", str(owl), "--queries", *queries, "-O", str(tmp_path)]
    ).returncode


# The full-file robot verify doesn't scale (male-cns's ~2 GB OWL exceeds heap), so
# the same checks run here on a small sample of EVERY generated file.
@pytest.mark.skipif(shutil.which("robot") is None, reason="ROBOT not on PATH")
def test_hemibrain_sample_passes_shape_checks(tmp_path, monkeypatch):
    _counts, _r2r, _report_dir, ofns = _hemibrain(monkeypatch, tmp_path)
    assert "ObjectPropertyAssertion(" in ofns[0].read_text()
    assert _shape_check(ofns[0], tmp_path) == 0


@pytest.mark.skipif(shutil.which("robot") is None, reason="ROBOT not on PATH")
def test_male_cns_sample_passes_shape_checks(tmp_path, monkeypatch):
    _report_dir, ofns = _male_cns(monkeypatch, tmp_path)
    for ofn in ofns:
        assert _shape_check(ofn, tmp_path) == 0


@pytest.mark.skipif(shutil.which("robot") is None, reason="ROBOT not on PATH")
def test_flywire_sample_passes_shape_checks(tmp_path, monkeypatch):
    # flywire pre edges carry `downstream` (no Tbars) -> the relaxed
    # pre-edge-missing-count check must still pass.
    _report_dir, ofns = _flywire(monkeypatch, tmp_path)
    assert _shape_check(ofns[0], tmp_path) == 0


def _totals_shape_check(ofn, tmp_path):
    owl = tmp_path / (ofn.stem + "_totals.owl")
    subprocess.run(["robot", "convert", "-i", str(ofn), "-f", "owl", "-o", str(owl)], check=True)
    return subprocess.run(
        ["robot", "verify", "-i", str(owl), "--queries",
         os.path.join(SPARQL, "synapse-total-violation.sparql"), "-O", str(tmp_path)]
    ).returncode


@pytest.mark.skipif(shutil.which("robot") is None, reason="ROBOT not on PATH")
def test_totals_sample_passes_shape_check(tmp_path, monkeypatch):
    _counts, _r2r, _report_dir, ofns = _hemibrain(monkeypatch, tmp_path)
    assert _totals_shape_check(ofns[0].parent / "synapse_totals_hemibrain.ofn", tmp_path) == 0


@pytest.mark.skipif(shutil.which("robot") is None, reason="ROBOT not on PATH")
def test_totals_shape_check_rejects_nonpositive_and_stray(tmp_path):
    for name, ax in (("zero", 'AnnotationAssertion(<%stotal_upstream> <%s> "0"^^xsd:integer)'),
                     ("stray", 'AnnotationAssertion(<%supstream> <%s> "5"^^xsd:integer)')):
        ofn = tmp_path / (name + ".ofn")
        ofn.write_text("Prefix(xsd:=<http://www.w3.org/2001/XMLSchema#>)\n" + _wrap(
            "Declaration(NamedIndividual(<%s>))" % _N, ax % (CUSTOM, _N)))
        assert _totals_shape_check(ofn, tmp_path) != 0


# --- the checks must REJECT the shapes we no longer emit -------------------
_N = "http://virtualflybrain.org/reports/VFB_00000001"   # a VFB neuron
_R = "http://virtualflybrain.org/reports/VFB_00000002"   # a VFB region individual
_FBBT = "http://purl.obolibrary.org/obo/FBbt_00007011"   # a class (not allowed as target)
_UP = CUSTOM + "upstream"


def _wrap(*axioms):
    return "Ontology(<http://x>\n" + "\n".join(axioms) + "\n)\n"


# each snippet is a shape we used to allow (or a bug) — verify must now flag it
_BAD_SHAPES = {
    # class form: ClassAssertion over a restriction (annotatedProperty=rdf:type,
    # target is anonymous) — dropped entirely
    "class_form": _wrap(
        "Declaration(NamedIndividual(<%s>))" % _N,
        "Declaration(Class(<%s>))" % _FBBT,
        "Declaration(ObjectProperty(<%s>))" % m.RO_POST,
        'ClassAssertion(Annotation(<%s> "5"^^xsd:float) '
        'ObjectSomeValuesFrom(<%s> <%s>) <%s>)' % (_UP, m.RO_POST, _FBBT, _N)),
    # a re-introduced per-side count bucket
    "side_bucket": _wrap(
        "Declaration(NamedIndividual(<%s>))" % _N,
        "Declaration(NamedIndividual(<%s>))" % _R,
        "Declaration(ObjectProperty(<%s>))" % m.RO_POST,
        'ObjectPropertyAssertion(Annotation(<%supstream_left> "5"^^xsd:float) '
        '<%s> <%s> <%s>)' % (CUSTOM, m.RO_POST, _N, _R)),
    # an edge to an FBbt class instead of a VFB region individual
    "nonvfb_target": _wrap(
        "Declaration(NamedIndividual(<%s>))" % _N,
        "Declaration(Class(<%s>))" % _FBBT,
        "Declaration(ObjectProperty(<%s>))" % m.RO_POST,
        'ObjectPropertyAssertion(Annotation(<%s> "5"^^xsd:float) '
        '<%s> <%s> <%s>)' % (_UP, m.RO_POST, _N, _FBBT)),
}


@pytest.mark.skipif(shutil.which("robot") is None, reason="ROBOT not on PATH")
@pytest.mark.parametrize("name", sorted(_BAD_SHAPES))
def test_shape_checks_reject_disallowed_forms(name, tmp_path):
    ofn = tmp_path / (name + ".ofn")
    ofn.write_text(_BAD_SHAPES[name])
    assert _shape_check(ofn, tmp_path) != 0

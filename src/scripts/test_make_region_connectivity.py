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
    "pre-edge-missing-tbars",
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
        assert entry["site"] and entry["roi_dataset"]


# --- ROI label parsing (regression for the backslash bug) ------------------
def test_roi_from_label_strips_suffix():
    assert m._roi_from_label("SIP(R) on JRC_FlyEM_Hemibrain") == "SIP(R)"


def test_roi_from_label_strips_backslash_escaped_apostrophe():
    # VFB stores a\'L(R); neuprint uses a'L(R)
    assert m._roi_from_label("a\\'L(R) on JRC_FlyEM_Hemibrain") == "a'L(R)"


def test_roi_from_label_handles_list():
    assert m._roi_from_label(["LH(R) on X"]) == "LH(R)"


# --- edge counts (primary gates the edge) ----------------------------------
POST = (m.RO_POST, "upstream", "upstream", ())
PRE = (m.RO_PRE, "pre", "Tbars", (("downstream", "downstream"),))


def test_edge_counts_postsynaptic():
    assert m._edge_counts({"upstream": 5}, *POST[1:]) == {"upstream": 5}
    assert m._edge_counts({"upstream": 0}, *POST[1:]) is None


def test_edge_counts_presynaptic_requires_tbars():
    # T-bars present -> edge; downstream included only when > 0
    assert m._edge_counts({"pre": 3, "downstream": 7}, *PRE[1:]) == {"Tbars": 3, "downstream": 7}
    assert m._edge_counts({"pre": 3, "downstream": 0}, *PRE[1:]) == {"Tbars": 3}
    # no T-bars -> no presynaptic edge even if downstream > 0 (the male-cns case)
    assert m._edge_counts({"pre": 0, "downstream": 4}, *PRE[1:]) is None


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


# --- integration: whole pipeline on saved samples (no live queries) ---------
def _run_main(dataset, tmp_path, monkeypatch, counts, b2n, np_dataset, roi_group, r2r):
    """Run main() fully offline: every network call is stubbed with saved data.
    Returns (report_dir, [produced OFN paths])."""
    monkeypatch.setattr(m, "get_neuprint_dataset", lambda site, url=None: np_dataset)
    monkeypatch.setattr(m, "get_bodyid_to_neuron", lambda site, url=None: b2n)
    monkeypatch.setattr(m, "get_roi_to_region", lambda ds, url=None: r2r)
    monkeypatch.setattr(m, "fetch_roi_counts", lambda client, ids: counts)
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

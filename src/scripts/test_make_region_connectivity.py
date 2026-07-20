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
    "nonvfb-neuron",
]


# --- helpers ---------------------------------------------------------------
def test_obo_iri_expands_curie():
    assert m._obo_iri("FBbt:00045032") == "http://purl.obolibrary.org/obo/FBbt_00045032"


def test_obo_iri_passthrough_iri():
    iri = "http://virtualflybrain.org/reports/VFB_00101436"
    assert m._obo_iri(iri) == iri


def test_lit_is_float_typed_integer_string():
    # counts are whole numbers but typed xsd:float to match pdb
    assert m._lit(23.0) == '"23"^^<http://www.w3.org/2001/XMLSchema#float>'
    assert m._lit(1665) == '"1665"^^<http://www.w3.org/2001/XMLSchema#float>'


# --- ROI label parsing (regression for the backslash bug) ------------------
def test_roi_from_label_strips_suffix():
    assert m._roi_from_label("SIP(R) on JRC_FlyEM_Hemibrain") == "SIP(R)"


def test_roi_from_label_strips_backslash_escaped_apostrophe():
    # VFB stores a\'L(R); neuprint uses a'L(R)
    assert m._roi_from_label("a\\'L(R) on JRC_FlyEM_Hemibrain") == "a'L(R)"


def test_roi_from_label_handles_list():
    assert m._roi_from_label(["LH(R) on X"]) == "LH(R)"


# --- annotation emission ---------------------------------------------------
def test_annotations_only_positive_counts():
    anns = m._annotations({"upstream": 5}, m.POST_ANNOTATIONS)
    assert anns == ['Annotation(<%supstream> "5"^^<http://www.w3.org/2001/XMLSchema#float>)' % CUSTOM]
    assert m._annotations({"upstream": 0}, m.POST_ANNOTATIONS) == []


def test_annotations_presynaptic_drops_zero_downstream():
    anns = m._annotations({"pre": 3, "downstream": 0}, m.PRE_ANNOTATIONS)
    joined = " ".join(anns)
    assert "custom/Tbars" in joined
    assert "downstream" not in joined


# --- axiom shape -----------------------------------------------------------
def test_relation_axiom_individual_is_object_property_assertion():
    ax = m._relation_axiom("n", "r", "individual", m.RO_POST, ['Annotation(x)'])
    assert ax.startswith("ObjectPropertyAssertion(")
    assert "<n>" in ax and "<r>" in ax


def test_relation_axiom_class_is_class_assertion():
    ax = m._relation_axiom("n", "FBbt", "class", m.RO_PRE, ['Annotation(x)'])
    assert ax.startswith("ClassAssertion(")
    assert "ObjectSomeValuesFrom(" in ax


# --- build_axioms ----------------------------------------------------------
def _counts():
    return pd.DataFrame([
        {"bodyId": 1, "roi": "SIP(R)", "upstream": 5, "pre": 2, "downstream": 10},
        {"bodyId": 1, "roi": "SCL(R)", "upstream": 3, "pre": 0, "downstream": 0},
        {"bodyId": 1, "roi": "SNP(R)", "upstream": 9, "pre": 9, "downstream": 9},  # unmapped
        {"bodyId": 99, "roi": "SIP(R)", "upstream": 1, "pre": 1, "downstream": 1},  # unknown neuron
    ])


_MAP = {
    "SIP(R)": ("http://vfb/sip", "individual"),
    "SCL(R)": ("http://vfb/scl", "individual"),
}


def test_build_axioms_counts_and_filtering():
    axioms, props = m.build_axioms(_counts(), {1: "http://vfb/n1"}, _MAP)
    # (1,SIP): post+pre = 2 ; (1,SCL): post only = 1 ; unmapped + unknown neuron dropped
    assert len(axioms) == 3
    assert props == {m.RO_POST, m.RO_PRE}


def test_build_axioms_scl_emits_only_postsynaptic():
    axioms, _ = m.build_axioms(_counts(), {1: "http://vfb/n1"}, _MAP)
    scl = [a for a in axioms if "http://vfb/scl" in a]
    assert len(scl) == 1
    assert m.RO_POST in scl[0] and m.RO_PRE not in scl[0]


# --- reports & ontology ----------------------------------------------------
def test_write_roi_report_flags_unmapped(tmp_path):
    out = tmp_path / "roi_mapping.tsv"
    m.write_roi_report(str(out), _counts(), _MAP)
    lines = out.read_text().splitlines()
    assert lines[0] == "roi\tkind\ttarget\tneurons"
    body = {ln.split("\t")[0]: ln.split("\t") for ln in lines[1:]}
    assert body["SNP(R)"][1] == "unmapped"
    assert body["SIP(R)"][1] == "individual"
    # unmapped ROIs float to the top
    assert lines[1].split("\t")[1] == "unmapped"


def test_write_ontology_header_has_dataset_and_version(tmp_path):
    out = tmp_path / "x.ofn"
    m.write_ontology(str(out), "hemibrain", "hemibrain:v1.2.1", ["ClassAssertion(a b)"], {m.RO_POST})
    text = out.read_text()
    assert 'dcterms:identifier "hemibrain"' in text
    assert 'dcterms:source "hemibrain:v1.2.1"' in text
    assert 'owl:versionInfo "v1.2.1"' in text
    assert "ClassAssertion(a b)" in text


# --- integration: whole pipeline on a saved sample (no live queries) --------
def _load_sample():
    counts = pd.read_csv(os.path.join(FIX, "neuprint_roi_counts_sample.tsv"), sep="\t")
    b2n = {int(k): v for k, v in json.load(open(os.path.join(FIX, "bodyid_to_neuron.json"))).items()}
    r2r = {k: tuple(v) for k, v in json.load(open(os.path.join(FIX, "roi_to_region.json"))).items()}
    return counts, b2n, r2r


def _run_pipeline(tmp_path, monkeypatch):
    """Run main() end-to-end against the saved sample; return (ofn, roi_report)."""
    counts, b2n, r2r = _load_sample()
    monkeypatch.setattr(m, "get_neuprint_dataset", lambda site: "hemibrain:v1.2.1")
    monkeypatch.setattr(m, "get_bodyid_to_neuron", lambda site: b2n)
    monkeypatch.setattr(m, "get_roi_to_region", lambda ds: r2r)
    monkeypatch.setattr(m, "fetch_roi_counts", lambda ds, tok, ids: counts)
    token = tmp_path / "token.txt"
    token.write_text("dummy")
    ofn = tmp_path / "region_connectivity_hemibrain.ofn"
    roi = tmp_path / "roi_mapping.tsv"
    monkeypatch.setattr(sys, "argv", ["prog", "hemibrain", str(token), str(ofn), str(roi)])
    m.main()
    return ofn, roi


def test_pipeline_from_sample(tmp_path, monkeypatch):
    counts, _, r2r = _load_sample()
    ofn, roi = _run_pipeline(tmp_path, monkeypatch)
    text = ofn.read_text()

    # header provenance
    assert 'dcterms:identifier "hemibrain"' in text
    assert 'owl:versionInfo "v1.2.1"' in text

    # edge counts match what the sample should yield (+1 for each Declaration)
    mapped = counts[counts["roi"].isin(r2r)]
    assert text.count("<%s>" % m.RO_POST) == int((mapped["upstream"] > 0).sum()) + 1
    assert text.count("<%s>" % m.RO_PRE) == int((mapped["pre"] > 0).sum()) + 1

    # unmapped neuprint super-ROIs are flagged and dropped (not emitted)
    kinds = {ln.split("\t")[0]: ln.split("\t")[1]
             for ln in roi.read_text().splitlines()[1:]}
    for sup in ("INP", "NotPrimary", "OL(R)"):
        assert kinds.get(sup) == "unmapped"


@pytest.mark.skipif(shutil.which("robot") is None, reason="ROBOT not on PATH")
def test_sample_owl_converts_and_passes_shape_checks(tmp_path, monkeypatch):
    ofn, _ = _run_pipeline(tmp_path, monkeypatch)
    owl = tmp_path / "sample.owl"
    subprocess.run(["robot", "convert", "-i", str(ofn), "-f", "owl", "-o", str(owl)], check=True)
    queries = [os.path.join(SPARQL, "%s-violation.sparql" % c) for c in RC_CHECKS]
    result = subprocess.run(
        ["robot", "verify", "-i", str(owl), "--queries", *queries, "-O", str(tmp_path)]
    )
    assert result.returncode == 0  # non-zero if any shape violation

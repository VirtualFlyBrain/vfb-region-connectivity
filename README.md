
![Build Status](https://github.com/VirtualFlyBrain/vfb-region-connectivity/actions/workflows/qc.yml/badge.svg)
# VFB Region Connectivity

Produces **region (neuropil) synaptic connectivity** OWL files for connectomics
datasets, for side-loading into [Virtual Fly Brain](https://virtualflybrain.org)
(VFB). One file is produced per dataset:

```
region_connectivity_<dataset>.owl
```

Each file records, for every neuron that has a VFB individual, how many synapses
that neuron has in each neuropil region — the per-region tallies that neuprint
reports in `roiInfo`.

## Representation

For each neuron × region, the counts are attached as **axiom annotations** in the
`n2o` `custom:` namespace (`http://n2o.neo/custom/`), so that
[`owl2neo4jcsv`](https://github.com/VirtualFlyBrain/neo4j2owl) turns each axiom
into a Neo4j edge whose properties are the counts. This matches the
representation already present in VFB's production graph (`pdb`):

| Relation | IRI | Direction | Count annotations |
|---|---|---|---|
| has presynaptic terminals in | `RO:0002113` | outputs | `Tbars` (T-bars / `pre`), `downstream` |
| has postsynaptic terminal in | `RO:0002110` | inputs | `upstream` |

Counts are typed `xsd:float` in the OWL (pdb stores them as float, e.g. `545.0`);
the reports cast them back to integers for readability. A postsynaptic
(`RO:0002110`) axiom is emitted when the neuron has any inputs in the region; a
presynaptic (`RO:0002113`) axiom when it has any T-bars there.

There is currently **no synapse-count threshold**: every region with at least one
synapse of the relevant polarity (count ≥ 1) is emitted. This matches the existing
VFB data, whose region counts also go down to 1. (If filtering is ever wanted, it
would be a `min_count` gate in `_annotations` in the generator.)

A region is emitted as **either**:

* an `ObjectPropertyAssertion` to a VFB region **individual** — for datasets that
  have per-dataset painted-domain individuals (e.g. hemibrain, whose regions are
  the `Xu2020roi` "painted domains", labelled `"<ROI> on JRC_FlyEM_Hemibrain"`); or
* a `ClassAssertion` over `ObjectSomeValuesFrom(...)` to an FBbt **class** — for
  datasets that have no region individuals.

Both forms convert to the same Neo4j edge shape.

The ontology header records the provenance: `dcterms:identifier` (the dataset
key, e.g. `hemibrain`), `owl:versionInfo` (the version, e.g. `v1.2.1`),
`dcterms:source` (the full neuprint dataset, e.g. `hemibrain:v1.2.1`), plus a
`dcterms:title` and `rdfs:comment`.

Example (functional syntax):

```
ObjectPropertyAssertion(
    Annotation(custom:upstream "245"^^xsd:float)
    obo:RO_0002110 VFB:VFB_jrchjz57 VFB:VFB_00101436)          # -> SIP(R)

ObjectPropertyAssertion(
    Annotation(custom:Tbars "159"^^xsd:float)
    Annotation(custom:downstream "1665"^^xsd:float)
    obo:RO_0002113 VFB:VFB_jrchjz57 VFB:VFB_00101436)
```

## Data source

Counts are pulled fresh from [neuprint](https://neuprint.janelia.org) with
`neuprint-python`. The neuprint dataset (e.g. `hemibrain:v1.2.1`) is **derived
from the VFB Site node** (`dataset=...` in its `link_base`), so the data always
matches the version VFB is pinned to rather than a hard-coded string.

The neuron and region identifiers come from `pdb`:

* **bodyId → VFB neuron IRI** — via the `database_cross_reference` accession to the
  dataset's Site node.
* **neuprint ROI → VFB region IRI** — region individuals `has_source` the dataset's
  ROI `DataSet`; the ROI name is the region label with the `" on <template>"`
  suffix removed (equal to the neuprint ROI name).

Only ROIs that resolve to a VFB region are emitted; this naturally drops
neuprint super-ROIs (e.g. `SNP(R)`, `INP`, `NotPrimary`), reproducing the region
set already in pdb.

## Building

Requires Docker (ODK image) and a neuprint token.

```bash
# put your neuprint token (single line) at the repo root
cp ~/neuprint_token.txt ./neuprint_token.txt

cd src/ontology
sh run.sh make                              # default goal: build all datasets in DATASETS
sh run.sh make ../../region_connectivity_hemibrain.owl   # build one dataset
sh run.sh make refresh                      # force a fresh neuprint pull for all datasets
```

Each build writes two files to the repo root: `region_connectivity_<dataset>.owl`
(the full side-load, ~160 MB for hemibrain) and `region_connectivity_<dataset>.owl.gz`
(~3 MB). The plain `.owl` is git-ignored (GitHub blocks files >100 MB); **the
`.owl.gz` is what is committed / released**, and the dumps pipeline decompresses
it. A dataset is only rebuilt when its `.owl` is missing — use `make refresh` (or
delete the file) to re-pull from neuprint.

## Reports & validation

```bash
cd src/ontology
sh run.sh make region_connectivity_reports   # term lists, summaries, shape checks
sh run.sh make obo_reports                    # optional ROBOT report (see note)
```

Written to `src/ontology/reports/` (per-dataset reports in a `<dataset>/`
subfolder; the merged term list at the root):

* **`<dataset>/terms.txt`** and a merged **`region_connectivity-terms.txt`** —
  sorted lists of the neurons and regions that appear. Diff these between
  versions to see which entities were added or removed.
* **`<dataset>/summary.tsv`** — the dataset/version/source (read from the OWL
  header) plus entity/edge counts and count totals (neurons, regions, pre/post
  edges, total upstream/Tbars/downstream). Diff to spot large shifts (e.g. a
  dropped ROI mapping or a neuprint re-export).
* **`<dataset>/roi_mapping.tsv`** — every neuprint ROI seen in the data, what it
  maps to (`kind` = `individual` VFB region / `class` FBbt term), the target IRI,
  and how many neurons touch it. ROIs with no VFB mapping are flagged `unmapped`
  (they are dropped from the OWL). Most unmapped ROIs are neuprint super-ROIs
  (`SNP(R)`, `INP`, `OL(R)`, `NotPrimary`, …) that would double-count; a large
  unmapped primary ROI signals a missing mapping. Written at build time (needs
  the neuprint pull), unlike the other reports.
* **`<dataset>/verify.txt`** — `PASS` marker; the build
  **fails** if any shape check in [`src/sparql`](src/sparql) returns violations:
  every postsynaptic (`RO:0002110`) edge has an `upstream`, every presynaptic
  (`RO:0002113`) edge has a `Tbars`, counts are positive, counts only sit on the
  two expected relations, and every edge source is a VFB neuron individual.

> The standard ODK ROBOT report (`<dataset>/obo-report.tsv`, via `make
> obo_reports`) is of limited use here: a side-load edge file intentionally
> carries no labels/definitions for external terms, so the report is almost
> entirely `missing_label`. `obo_reports` filters those rows out.

## Tests & CI

`src/scripts/test_make_region_connectivity.py` covers, with **no live queries**:

* **unit tests** of the pure logic — ROI-label parsing (incl. the backslash fix),
  annotation/axiom emission, edge filtering, report generation, header metadata;
* an **offline pipeline integration test** — the network calls are stubbed with a
  *saved* neuprint/pdb sample in [`src/scripts/test_data/`](src/scripts/test_data)
  (a few neurons' ROI counts + the bodyId→neuron and ROI→region maps), and the
  real transform runs end-to-end through `main()`; a second test ROBOT-converts
  that output and asserts it passes the SPARQL shape checks.

```bash
pip install -r src/scripts/requirements.txt pytest && pytest src/scripts
```

The [`qc.yml`](.github/workflows/qc.yml) GitHub Actions workflow runs two jobs on
every push/PR (both in the ODK image, no neuprint token needed):

* **tests** — the pytest suite above (unit + offline pipeline integration).
* **validate-artifacts** — gunzips each committed `region_connectivity_*.owl.gz`
  and runs the SPARQL shape checks against it (via the Makefile `verify` target),
  so a malformed committed artifact fails CI.

To refresh the sample after a neuprint/schema change, regenerate the three files
in `test_data/` from a small fixed set of bodyIds.

## Adding a dataset

1. Add a `DATASETS` entry in
   [`src/scripts/make_region_connectivity.py`](src/scripts/make_region_connectivity.py)
   with the VFB `site` short_form and either a `roi_dataset` (region individuals)
   or a `roi_mapping_tsv` (ROI → FBbt class curation).
2. Add the dataset key to `DATASETS` in
   [`src/ontology/VFB_region_connectivity.Makefile`](src/ontology/VFB_region_connectivity.Makefile).

## Side-loading

The output files are consumed by
[`vfb-pipeline-dumps`](https://github.com/VirtualFlyBrain/vfb-pipeline-dumps) as
edge-only side-loads (kept out of the reasoned/merged graph). They are named
`region_connectivity_*` to distinguish them from the neuron-to-neuron
`connectome_*` side-loads; the dumps pipeline's `PDB_EXTERNAL_ONTS` list needs to
include `region_connectivity_*.owl`.

## Repository layout

| Path | What it is |
|---|---|
| `region_connectivity_<dataset>.owl.gz` | The committed side-load artifact(s) (repo root). The uncompressed `.owl` is git-ignored. |
| `src/scripts/make_region_connectivity.py` | The generator: neuprint pull → mappings → OWL + ROI report. |
| `src/scripts/cypher_query.py` | Helper for querying VFB's Neo4j (`pdb`). |
| `src/scripts/requirements.txt` | Python dependencies. |
| `src/scripts/test_make_region_connectivity.py` | Unit + offline pipeline tests. |
| `src/scripts/test_data/` | Saved neuprint/pdb sample used by the pipeline test. |
| `src/ontology/VFB_region_connectivity.Makefile` | Custom build & report targets (edit this, not `Makefile`). |
| `src/ontology/Makefile`, `run.sh` | ODK-managed wrappers (do not edit). |
| `src/ontology/reports/` | Generated reports — per-dataset in `<dataset>/`, merged term list at the root. |
| `src/sparql/region_connectivity_*.sparql` | Term-list and summary report queries. |
| `src/sparql/*-violation.sparql` | Shape-validation checks (the repo's own five, plus ODK defaults). |
| `.github/workflows/qc.yml` | CI: tests + committed-artifact validation. |

## Acknowledgements

This repository was created using the
[Ontology Development Kit (ODK)](https://github.com/INCATools/ontology-development-kit).

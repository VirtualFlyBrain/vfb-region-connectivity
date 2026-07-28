
![Build Status](https://github.com/VirtualFlyBrain/vfb-region-connectivity/actions/workflows/qc.yml/badge.svg)
# VFB Region Connectivity

Produces one **region (neuropil) synaptic connectivity** OWL file per
connectomics dataset, for side-loading into
[Virtual Fly Brain](https://virtualflybrain.org) (VFB):

```
region_connectivity_<dataset>.owl
```

Each file records, for every neuron that has a VFB individual, how many synapses
that neuron has in each neuropil region — the per-region tallies that neuprint
reports in `roiInfo`. Even a whole-CNS dataset (male-cns, banc) compresses to a
few tens of MB as a single file, so no splitting is needed.

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

Each edge is an `ObjectPropertyAssertion` from the neuron to a VFB region
**individual**, carrying the count total(s). The individual is side-specific, so
**laterality is carried by the target** — no per-side annotation is stored; it is
recovered downstream by joining the region's side (region individuals carry it
via `part_of` in the KB) with the neuron's soma side (VFB `has_soma_location`).

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

The neuron and region identifiers come from VFB's Neo4j. By default this is the
**knowledge base (`kb`)**, not `pdb`: dataset updates land in the kb first, so a
new dataset version can be built before it reaches `pdb`. (An `endpoint: pdb`
per-dataset override is available.) Neurons are the non-deprecated `Individual`s
cross-referenced to the Site — i.e. the cell-typed subset VFB loads.

* **bodyId → VFB neuron IRI** — via the `database_cross_reference` accession to the
  dataset's Site node.
* **neuprint ROI → VFB region individual** (`roi_dataset`): the per-side
  individuals that `has_source` the dataset's ROI `DataSet` (hemibrain
  `Xu2020roi`, male-cns `Berg2025aroi`, banc `Bates2026roi`). The ROI name is the
  region label minus the `" on <template>"` suffix.

Only ROIs that resolve to a region individual are emitted; this drops ROIs with no mapping
(neuprint super-ROIs like `SNP(R)`/`INP`, and male-cns's optic-lobe column
sub-ROIs — whose synapses are already counted at their parent region), which are
listed in the ROI-mapping report.

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
sh run.sh make region_connectivity_reports   # term lists + summaries
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
* **`<dataset>/roi_mapping.tsv`** — every neuprint ROI seen in the data, its
  `group`, the VFB region-individual `target` IRI (blank if unmapped), and how
  many neurons touch it. The group is `unmapped` for any ROI with no region
  individual (otherwise its brain/vnc group). Unmapped ROIs are dropped from the
  OWL — mostly neuprint super-ROIs (`SNP(R)`, `INP`, …) that would double-count,
  or optic-lobe columns; a large unmapped primary ROI signals a missing region
  individual. Written at build time (needs the neuprint pull), unlike the other
  reports.
**Shape validation** ([`src/sparql/*-violation.sparql`](src/sparql)) locks the
output to exactly the shape we now emit — a direct `ObjectPropertyAssertion` from
a VFB neuron to a VFB region individual, annotated with the three totals:

* every postsynaptic (`RO:0002110`) edge has an `upstream`, every presynaptic
  (`RO:0002113`) edge has a `Tbars`, and counts are positive;
* counts sit **only** on `upstream`/`Tbars`/`downstream` (no stray or per-side
  annotation) and **only** on those two relations (any `rdf:type`/class-form
  axiom is rejected);
* every edge **source** is a VFB neuron and every edge **target** is a VFB region
  individual (an FBbt class or anonymous restriction is rejected).

This is **not** run over the full file (`robot verify` doesn't fit in memory at
connectome scale); the same checks run on small samples in the test suite (below,
which also asserts the disallowed shapes are rejected), and the invariants are
guaranteed by construction in `build_axioms`.

> The standard ODK ROBOT report (`<dataset>/obo-report.tsv`, via `make
> obo_reports`) is of limited use here: a side-load edge file intentionally
> carries no labels/definitions for external terms, so the report is almost
> entirely `missing_label`. `obo_reports` filters those rows out.

## Tests & CI

`src/scripts/test_make_region_connectivity.py` covers, with **no live queries**:

* **unit tests** of the pure logic — ROI-label parsing (incl. the backslash fix),
  edge-count gating, aggregation, report generation, header metadata;
* **offline pipeline tests** — network calls are stubbed with *saved* samples in
  [`src/scripts/test_data/`](src/scripts/test_data), and the real transform runs
  end-to-end through `main()` for hemibrain and male-cns; each output is
  ROBOT-converted and run through the `src/sparql/*-violation.sparql` shape checks.

```bash
pip install -r src/scripts/requirements.txt pytest && pytest src/scripts
```

The [`qc.yml`](.github/workflows/qc.yml) GitHub Actions workflow runs the pytest
suite above on every push/PR (in the ODK image, no neuprint token needed). The
full committed artifacts are not shape-checked in CI — that doesn't fit in a
runner's memory at connectome scale; the sample checks cover both forms.

To refresh the sample after a neuprint/schema change, regenerate the three files
in `test_data/` from a small fixed set of bodyIds.

## Adding a dataset

Add one entry to
[`src/ontology/dataset_config.yaml`](src/ontology/dataset_config.yaml) — the
single source of truth, read by both the generator and the Makefile:

```yaml
<dataset-key>:
  site: <VFB Site short_form>          # the connectome data source
  roi_dataset: <VFB DataSet short_form>  # the region individuals
  # endpoint: pdb                      # optional; defaults to kb
```

The top-level key is the `<dataset>` build arg and the output-file suffix. No
code or Makefile change is needed — `make region_connectivity` picks it up.

The region individuals must exist in VFB first (loaded from the connectome's
painted-domain ROIs). ROIs with no individual are dropped and flagged in
`roi_mapping.tsv`.

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
| `src/scripts/cypher_query.py` | Helper for querying VFB's Neo4j (`kb`/`pdb`). |
| `src/scripts/requirements.txt` | Python dependencies. |
| `src/scripts/test_make_region_connectivity.py` | Unit + offline pipeline tests. |
| `src/scripts/test_data/` | Saved neuprint/pdb sample used by the pipeline test. |
| `src/ontology/dataset_config.yaml` | Per-dataset config (site / roi_dataset / endpoint) — read by the generator and the Makefile. |
| `src/ontology/VFB_region_connectivity.Makefile` | Custom build & report targets (edit this, not `Makefile`). |
| `src/ontology/Makefile`, `run.sh` | ODK-managed wrappers (do not edit). |
| `src/ontology/reports/` | Generated reports — per-dataset in `<dataset>/`, merged term list at the root. |
| `src/sparql/region_connectivity_*.sparql` | Term-list and summary report queries. |
| `src/sparql/*-violation.sparql` | Shape-validation checks that lock the output to the emitted individual form. |
| `.github/workflows/qc.yml` | CI: tests + committed-artifact validation. |

## Acknowledgements

This repository was created using the
[Ontology Development Kit (ODK)](https://github.com/INCATools/ontology-development-kit).

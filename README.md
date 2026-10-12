
![Build Status](https://github.com/VirtualFlyBrain/vfb-region-connectivity/actions/workflows/qc.yml/badge.svg)
# VFB Region Connectivity

Produces one **region (neuropil) synaptic connectivity** OWL file per
connectomics dataset, for side-loading into
[Virtual Fly Brain](https://virtualflybrain.org) (VFB):

```
region_connectivity_<dataset>.owl
```

Each file records, for every neuron that has a VFB individual, how many synapses
that neuron has in each neuropil region. Counts come from the dataset's source —
either [neuprint](https://neuprint.janelia.org) (hemibrain, male-cns, manc,
optic-lobe) or a [flywire codex](https://codex.flywire.ai) CSV (flywire/FAFB).
Even a whole-CNS dataset compresses to a few tens of MB as a single file, so no
splitting is needed.

The same pull also produces one **whole-neuron synapse totals** file per dataset:

```
synapse_totals_<dataset>.owl
```

holding each neuron's total input and output synapse counts (see
[Whole-neuron synapse totals](#whole-neuron-synapse-totals)).

## Representation

For each neuron × region, the counts are attached as **axiom annotations** in the
`n2o` `custom:` namespace (`http://n2o.neo/custom/`), so that
[`owl2neo4jcsv`](https://github.com/VirtualFlyBrain/neo4j2owl) turns each axiom
into a Neo4j edge whose properties are the counts. This matches the
representation already present in VFB's production graph (`pdb`):

| Relation | IRI | Direction | Count annotations |
|---|---|---|---|
| has presynaptic terminals in | `RO:0002113` | outputs | `Tbars` (T-bars / `pre`) + `downstream` (neuprint); `downstream` only (flywire) |
| has postsynaptic terminal in | `RO:0002110` | inputs | `upstream` |

flywire has no T-bar annotation, so its presynaptic edge carries only
`downstream` (the output-synapse count); its postsynaptic edge carries `upstream`
(the input-synapse count).

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

## Whole-neuron synapse totals

`synapse_totals_<dataset>.owl` gives each VFB neuron its total synapse counts as
annotation assertions in the `custom:` namespace, which `owl2neo4jcsv` turns into
**node properties** (the same mechanism as the scRNAseq `cell_count`):

| Node property | Meaning | neuprint | flywire | banc | catmaid |
|---|---|---|---|---|---|
| `total_upstream` | input synapses | `:Neuron.upstream` (`post` if absent) | `per_neuron_neuropil_count_post_783`, summed | parquet rows where the neuron is post, summed | postsynaptic links |
| `total_downstream` | output synapses | `:Neuron.downstream` | `per_neuron_neuropil_count_pre_783`, summed | parquet rows where the neuron is pre, summed | — |
| `total_Tbars` | presynaptic sites | `:Neuron.pre` | — (no T-bars) | — | presynaptic connectors |

BANC and CATMAID are **totals-only** datasets (`totals_only: true`): they get a
`synapse_totals_<dataset>.owl` but no region file yet. BANC (v888, Site
`BANC888`) has region individuals (`Bates2026roi`) and per-neuropil counts, but
the parquet's neuropil names (`CB_SLP_R`, …) don't match the region labels yet.

* **banc** — `totals_parquet_url`: the v888
  `synapses_v1_human_readable_sizethresh3_connectioncountsperneuropil.parquet`,
  i.e. the v1 synapse set (size > 3) with no connection-count threshold — the
  same synapses VFB's BANC connectivity uses before its count >= 3 cut. Partners
  that are fragments or unsegmented (root 0) still count. (`banc_888_metrics`
  also has per-neuron totals, but from the v2 synapse set, so they wouldn't match
  the edges.)
* **catmaid** — the VFB-hosted CATMAID server and project are read from the
  Site's `link_base`; each skeleton's `compact-detail` connectors give its input
  links and presynaptic connectors, traced partner or not (one request per
  skeleton, 8 at a time). Synapses *made* would need every connector's partners,
  so there is no `total_downstream`. A live server has no version, so the build
  date is recorded as `dcterms:source` `<site>:<date>`.

```
AnnotationAssertion(custom:total_upstream VFB:VFB_jrchjz57 "1531"^^xsd:integer)   # illustrative value
```

These are the denominators for input-proportion connectivity weights (the share
of a neuron's input that comes from a given partner), planned for the VFB circuit
browser and class-level connectivity. They **cannot** be derived from the region edges or from VFB's
`synapsed_to` edges:

* mapped neuprint ROIs can overlap (hemibrain maps super-ROIs as well as their
  parts, e.g. `CRE(R)` and `CRE(-ROB,-RUB)(R)`), and unmapped ones are dropped,
  so region counts don't sum to a neuron total;
* flywire's `neuropil_synapse_table` only counts connections of **>= 5
  synapses**, so its region counts are thresholded;
* VFB's `synapsed_to` edges only connect neurons VFB has loaded and drop
  connections of 1 (or, for BANC, 2) synapses.

Every source counts **every synapse on the neuron, whatever the partner**,
unproofread fragments included. flywire's come from the FlyWire 783
per-neuron neuropil synapse counts on
[Zenodo 10676866](https://zenodo.org/records/10676866) (`totals_post_url` /
`totals_pre_url` in `dataset_config.yaml`; cleft score >= 50, the same synapse
set as VFB's flywire connectivity), not from codex, whose tables only cover
connections between proofread neurons. Each file's `rdfs:comment` records its
source. During the build, `report_totals` prints
how the totals compare with the mapped region input counts.

**Loading:** unlike the region files, these go into the main graph, not the
edge-only side-load. `vfb-pipeline-dumps` converts side-loads with
`owl2neo4jcsv … only_edges`, which would drop node properties. So they need
listing in `vfb-pipeline-collectdata`'s `vfb_fullontologies.txt` (e.g.
`http://virtualflybrain.org/data/VFB/OWL/synapse_totals_*`), like the scRNAseq
metadata files, and must **not** match `PDB_EXTERNAL_ONTS`.

## Data source

The count source is set per dataset (`source` in `dataset_config.yaml`):

* **`neuprint`** (default) — counts pulled fresh with `neuprint-python`. The
  neuprint dataset (e.g. `hemibrain:v1.2.1`) is **derived from the VFB Site node**
  (`dataset=...` in its `link_base`), so it always matches the version VFB is
  pinned to. Uses neuprint's `roiInfo` (`pre`/`upstream`/`downstream`).
* **`flywire`** — counts read from a flywire codex CSV (`csv_url`, the
  `neuropil_synapse_table` for release 783). One value per direction per region
  (`input synapses in <R>` → `upstream`, `output synapses in <R>` → `downstream`);
  no T-bars. The version (`flywire:783`) is derived from the Site's `data_version`.

The neuron and region identifiers always come from VFB's Neo4j. By default this
is the **knowledge base (`kb`)**, not `pdb`: dataset updates land in the kb first,
so a new dataset version can be built before it reaches `pdb`. (An `endpoint: pdb`
per-dataset override is available.) Neurons are the non-deprecated `Individual`s
cross-referenced to the Site — i.e. the cell-typed subset VFB loads.

* **bodyId / root_id → VFB neuron IRI** — via the `database_cross_reference`
  accession to the dataset's Site node.
* **region name → VFB region individual** (`roi_dataset`): the per-side
  individuals that `has_source` the dataset's ROI `DataSet` (hemibrain
  `Xu2020roi`, male-cns `Berg2025aroi`, manc `Takemura2023roi`, optic-lobe
  `Nern2024roi`, flywire `Dorkenwald2023roi`). The name is the region label minus
  the `" on <template>"` suffix.

Only regions that resolve to an individual are emitted; this drops names with no
mapping (neuprint super-ROIs like `SNP(R)`/`INP`, optic-lobe column sub-ROIs whose
synapses are already counted at their parent region, and flywire's `UNASGD`),
which are listed in the ROI-mapping report.

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
  (`RO:0002113`) edge has a presynaptic count (`Tbars` for neuprint, or
  `downstream` for flywire), and counts are positive;
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
* **offline pipeline tests** — network calls (and the flywire CSV read) are
  stubbed with *saved* samples in
  [`src/scripts/test_data/`](src/scripts/test_data), and the real transform runs
  end-to-end through `main()` for a neuprint dataset (hemibrain, male-cns) and the
  flywire CSV path; each output is ROBOT-converted and run through the
  `src/sparql/*-violation.sparql` shape checks (which also assert the disallowed
  shapes are rejected).

```bash
pip install -r src/scripts/requirements.txt pytest && pytest src/scripts
```

The [`qc.yml`](.github/workflows/qc.yml) GitHub Actions workflow runs the pytest
suite above on every push/PR (in the ODK image, no neuprint token needed). The
full committed artifacts are not shape-checked in CI — that doesn't fit in a
runner's memory at connectome scale; the sample checks cover every source path.

To refresh the sample after a neuprint/schema change, regenerate the three files
in `test_data/` from a small fixed set of bodyIds.

## Adding a dataset

Add one entry to
[`src/ontology/dataset_config.yaml`](src/ontology/dataset_config.yaml) — the
single source of truth, read by both the generator and the Makefile:

```yaml
<dataset-key>:
  site: <VFB Site short_form>            # the connectome data source
  roi_dataset: <VFB DataSet short_form>  # the region individuals
  # endpoint: pdb                        # optional; defaults to kb
  # source: flywire                      # optional; defaults to neuprint
  # csv_url: https://…/table.csv.gz      # required when source: flywire
  # totals_only: true                    # totals file only (no roi_dataset needed)
```

The top-level key is the `<dataset>` build arg and the output-file suffix. No
code or Makefile change is needed — `make region_connectivity` picks it up. A
neuprint dataset needs only `site` + `roi_dataset`; a flywire dataset also needs
`source: flywire` and the `csv_url` of its `neuropil_synapse_table` CSV (plus
`totals_post_url` / `totals_pre_url`). A `totals_only` dataset (source `banc` or
`catmaid`) needs no `roi_dataset`; the Makefile builds it with its own rule.

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
| `synapse_totals_<dataset>.owl.gz` | Whole-neuron synapse totals (node properties), loaded into the main graph. The uncompressed `.owl` is git-ignored. |
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
| `src/sparql/*-violation.sparql` | Shape-validation checks that lock the output to the emitted individual form (`synapse-total-violation.sparql` for the totals files). |
| `.github/workflows/qc.yml` | CI: tests + committed-artifact validation. |

## Acknowledgements

This repository was created using the
[Ontology Development Kit (ODK)](https://github.com/INCATools/ontology-development-kit).

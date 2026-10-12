
![Build Status](https://github.com/VirtualFlyBrain/vfb-region-connectivity/actions/workflows/qc.yml/badge.svg)
# VFB Region Connectivity

Produces, for each connectomics dataset in
[Virtual Fly Brain](https://virtualflybrain.org) (VFB), up to two OWL files:

```
region_connectivity_<dataset>.owl   # per neuron x neuropil region: synapse counts (edges)
synapse_totals_<dataset>.owl        # per neuron: whole-neuron synapse totals (node properties)
```

The region file records, for every neuron that has a VFB individual, how many
synapses it has in each neuropil region; it is side-loaded into VFB's graph as
edges. The totals file records each neuron's total input and output synapses —
the denominators for input-proportion connectivity weights — and is loaded into
the main graph as node properties (see [Loading](#loading)).

Every count includes **every synapse on the neuron, whatever the partner**
(unproofread fragments included) and has **no connection-count threshold**, so
counts mean the same thing across datasets.

| Dataset | VFB Site | Source | Region file | Totals file |
|---|---|---|---|---|
| `hemibrain` | `neuprint_JRC_Hemibrain_1point2point1` | neuprint `hemibrain:v1.2.1` | ✓ (`Xu2020roi`) | ✓ |
| `male-cns` | `male-cns_v1_0` | neuprint `male-cns:v1.0` | ✓ (`Berg2025aroi`) | ✓ |
| `manc` | `neuprint_JRC_Manc_1_2_1` | neuprint `manc:v1.2.1` | ✓ (`Takemura2023roi`) | ✓ |
| `optic-lobe` | `neuprint_JRC_OpticLobe_v1_0_1` | neuprint `optic-lobe:v1.0.1` | ✓ (`Nern2024roi`) | ✓ |
| `flywire` | `flywire783` | FlyWire 783 per-neuron neuropil counts ([Zenodo 10676866](https://zenodo.org/records/10676866)) | ✓ (`Dorkenwald2023roi`) | ✓ |
| `banc` | `BANC888` | BANC v888 v1 per-neuropil connection counts (BANC GCS bucket) | ✓ (`Bates2026roi`) | ✓ |
| `catmaid-l1em` | `catmaid_l1em` | VFB-hosted CATMAID server | — | ✓ |
| `catmaid-fafb` | `catmaid_fafb` | VFB-hosted CATMAID server | — | ✓ |

All datasets are configured in
[`src/ontology/dataset_config.yaml`](src/ontology/dataset_config.yaml). Even a
whole-CNS region file compresses to a few tens of MB, so no splitting is needed.

## Representation

### Region file

For each neuron × region, the counts are attached as **axiom annotations** in the
`n2o` `custom:` namespace (`http://n2o.neo/custom/`), so that
[`owl2neo4jcsv`](https://github.com/VirtualFlyBrain/neo4j2owl) turns each axiom
into a Neo4j edge whose properties are the counts. This matches the
representation already present in VFB's production graph (`pdb`):

| Relation | IRI | Direction | Count annotations |
|---|---|---|---|
| has presynaptic terminals in | `RO:0002113` | outputs | `Tbars` (T-bars / `pre`) + `downstream` (neuprint); `downstream` only (flywire, banc) |
| has postsynaptic terminal in | `RO:0002110` | inputs | `upstream` |

flywire and banc have no T-bar annotation, so their presynaptic edge carries only
`downstream` (the output-synapse count) and is gated on it; for neuprint it is
gated on T-bars. A postsynaptic (`RO:0002110`) axiom is emitted when the neuron
has any inputs in the region.

Counts are typed `xsd:float` in the OWL (pdb stores them as float, e.g. `545.0`);
the reports cast them back to integers for readability. There is **no
synapse-count threshold**: every region with at least one synapse of the
relevant polarity (count ≥ 1) is emitted. (If filtering is ever wanted, it would
be a gate in `_edge_counts` in the generator.)

Each edge is an `ObjectPropertyAssertion` from the neuron to a VFB region
**individual**, carrying the count total(s). The individual is side-specific, so
**laterality is carried by the target** — no per-side annotation is stored; it is
recovered downstream by joining the region's side (region individuals carry it
via `part_of` in the KB) with the neuron's soma side (VFB `has_soma_location`).

```
ObjectPropertyAssertion(
    Annotation(custom:upstream "245"^^xsd:float)
    obo:RO_0002110 VFB:VFB_jrchjz57 VFB:VFB_00101436)          # -> SIP(R)

ObjectPropertyAssertion(
    Annotation(custom:Tbars "159"^^xsd:float)
    Annotation(custom:downstream "1665"^^xsd:float)
    obo:RO_0002113 VFB:VFB_jrchjz57 VFB:VFB_00101436)
```

Mapped regions can overlap (hemibrain and male-cns include super-ROIs as well as
their parts, e.g. `CRE(R)` and `CRE(-ROB,-RUB)(R)`), so region counts must not be
summed into a neuron total — use the totals file.

### Totals file

Each VFB neuron gets its whole-neuron synapse totals as annotation assertions in
the `custom:` namespace, which `owl2neo4jcsv` turns into **node properties** (the
same mechanism as the scRNAseq `cell_count`). They are typed `xsd:integer`, and
zero totals are omitted.

| Node property | Meaning | neuprint | flywire | banc | catmaid |
|---|---|---|---|---|---|
| `total_upstream` | input synapses | `:Neuron.upstream` (`post` if absent) | post file, summed over neuropils | parquet rows where the neuron is post, summed | postsynaptic links |
| `total_downstream` | output synapses | `:Neuron.downstream` | pre file, summed over neuropils | parquet rows where the neuron is pre, summed | — |
| `total_Tbars` | presynaptic sites | `:Neuron.pre` | — (no T-bars) | — | presynaptic connectors |

```
AnnotationAssertion(custom:total_upstream VFB:VFB_jrchjz57 "1531"^^xsd:integer)   # illustrative value
```

These are the denominators for input-proportion connectivity weights (the share
of a neuron's input that comes from a given partner), for the VFB circuit
browser and class-level connectivity. They **cannot** be derived from the region
edges (regions overlap or are unmapped) or from VFB's `synapsed_to` edges, which
only connect neurons VFB has loaded and drop connections of 1 (or, for BANC, 2)
synapses.

### Provenance

Both files' ontology headers record `dcterms:identifier` (the dataset key, e.g.
`hemibrain`), `dcterms:source` (`<source>:<version>`, e.g. `hemibrain:v1.2.1`,
`flywire:783`, `banc:888`, or `<site>:<build date>` for CATMAID, whose live
servers have no version), `owl:versionInfo`, a `dcterms:title` and an
`rdfs:comment` that, for the totals file, says exactly what the counts include.

## Data sources

The source is set per dataset (`source` in `dataset_config.yaml`):

* **`neuprint`** (default) — counts pulled fresh with `neuprint-python`. The
  neuprint dataset (e.g. `hemibrain:v1.2.1`) is **derived from the VFB Site node**
  (`dataset=...` in its `link_base`), so it always matches the version VFB is
  pinned to. Region counts come from neuprint's per-ROI counts
  (`pre`/`upstream`/`downstream`), totals from the `:Neuron` node's own counts.
* **`flywire`** — the FlyWire 783 per-neuron neuropil synapse counts on
  [Zenodo 10676866](https://zenodo.org/records/10676866) (`counts_post_url` /
  `counts_pre_url`): synapses received (`upstream`) and made (`downstream`) per
  neuron and neuropil, cleft score >= 50 (the synapse set VFB's flywire
  connectivity uses). Synapses in no neuropil (`None`) are relabelled `UNASGD`,
  codex's name. These replaced codex's `neuropil_synapse_table`, which only
  counts connections of >= 5 synapses between proofread neurons. Zenodo times
  out about half of requests, so downloads retry (`_local_copy`). The version
  (`flywire:783`) is derived from the Site's `data_version`.
* **`banc`** — BANC's v888 per-neuropil connection counts parquet
  (`counts_parquet_url`,
  `synapses_v1_human_readable_sizethresh3_connectioncountsperneuropil.parquet`):
  the v1 synapse set (size > 3), the one VFB's BANC connectivity is built from,
  with no connection-count threshold; partners that are fragments or
  unsegmented (root 0) still count. Read one row group at a time (~177M rows).
  The version (`banc:888`) is derived from the Site's `data_version`. Its
  neuropil labels (`{area}_{neuropil}_{side}`, e.g. `CB_SLP_R`,
  `VNC_T1-ProNm_L`) are BANC's own renamings of the Ito (brain) and Court (VNC)
  meshes in its `region_outlines` layer — the meshes `Bates2026roi` was curated
  from — made by `synapses/v1.1/finalize_neuropil_assignments.py` in the BANC
  bucket. `banc_roi_name` inverts them: `T1-ProNm` → `ProNM-T1`, `ANm` →
  `ABDNM`, `AMNp` → `AMNP`, `MB-CA` → `MB_CA`, Court's position-split
  `T1/T2/T3-mVAC` → `mVAC`, and `IntTct`/`LTct` → the `(COURT)` versions. The
  sided region is used where VFB has one (`SLP_R`), else the unsided one (`GNG`,
  `ProNM-T1`), so both sides' synapses are summed on it. The partition has no
  glomeruli, tracts or nerves, and its VNC uses the Court meshes, so the MANC-atlas
  regions in `Bates2026roi` get no counts. (`banc_888_metrics` also has per-neuron
  totals, but from the v2 synapse set, so they wouldn't match VFB's edges.)
* **`catmaid`** (totals only) — the VFB-hosted CATMAID server and project are
  read from the Site's `link_base`; each skeleton's `compact-detail` connectors
  give its input links and the connectors it is presynaptic at, traced partner or
  not (one request per skeleton, 8 at a time). Synapses *made* would need every
  connector's partners fetched, so there is no `total_downstream`.

The neuron and region identifiers always come from VFB's Neo4j. By default this
is the **knowledge base (`kb`)**, not `pdb`: dataset updates land in the kb first,
so a new dataset version can be built before it reaches `pdb`. (An `endpoint: pdb`
per-dataset override is available.) Neurons are the non-deprecated `Individual`s
cross-referenced to the Site — i.e. the cell-typed subset VFB loads.

* **bodyId / root_id / skeleton id → VFB neuron IRI** — via the
  `database_cross_reference` accession to the dataset's Site node.
* **region name → VFB region individual** (`roi_dataset`): the per-side
  individuals that `has_source` the dataset's ROI `DataSet`. A region answers to
  its label minus the `"<Dataset>:"` prefix and `" on <template>"` suffix
  (`MANC:LegNp(T3)(R) on JRC2018UnisexVNC` → `LegNp(T3)(R)`), and to its synonyms
  (`synonyms`, or hemibrain's `has_exact_synonym`), which newer ROI datasets set
  to the exact source ROI name (BANC `LTct (COURT)` has synonym `LTct`). A name
  matching more than one region is ambiguous and dropped with a warning (BANC's
  Court and MANC tracts share synonyms; their labels still tell them apart).

Only regions that resolve to an individual are emitted; this drops names with no
mapping (neuprint `NotPrimary` and `*-unspecified` catch-alls, optic-lobe column
sub-ROIs whose synapses are already counted at their parent region, flywire's
`UNASGD`, BANC's `*_unassigned_*`), which are listed in the ROI-mapping report.

## Building

Requires Docker (ODK image), access to VFB's kb, and — for the neuprint datasets —
a neuprint token.

```bash
# put your neuprint token (single line) at the repo root
cp ~/neuprint_token.txt ./neuprint_token.txt

cd src/ontology
sh run.sh make                                         # build every dataset
sh run.sh make ../../region_connectivity_hemibrain.owl # one dataset (region + totals)
sh run.sh make ../../synapse_totals_catmaid-l1em.owl   # one totals-only dataset
sh run.sh make refresh                                 # force a fresh pull for all datasets
```

A region dataset's build writes `region_connectivity_<dataset>.owl` and
`synapse_totals_<dataset>.owl` to the repo root, both from one pull, plus a
`.owl.gz` of each; a `totals_only` dataset writes just the totals file. The plain
`.owl` files are git-ignored (GitHub blocks files >100 MB); **the `.owl.gz`
files are what is committed / released**. Each build also shape-checks a sample
of its output (see below). A dataset is only rebuilt when its `.owl` is missing —
use `make refresh` (or delete the file) to re-pull.

## Reports & validation

```bash
cd src/ontology
sh run.sh make region_connectivity_reports   # term lists + summaries
sh run.sh make obo_reports                    # optional ROBOT report (see note)
```

Written to `src/ontology/reports/` (per-dataset reports in a `<dataset>/`
subfolder; the merged term list at the root), and committed so versions can be
diffed in git:

* **`<dataset>/terms.txt`** and a merged **`region_connectivity-terms.txt`** —
  sorted lists of the neurons and regions that appear in the region files.
  Diff these between versions to see which entities were added or removed.
* **`<dataset>/summary.tsv`** — the region file's dataset/version/source (read
  from the OWL header) plus entity/edge counts and count totals (neurons,
  regions, pre/post edges, total upstream/Tbars/downstream). Diff to spot large
  shifts (e.g. a dropped ROI mapping or a source re-export).
* **`<dataset>/totals_summary.tsv`** — the same for the totals file: neurons with
  totals and the summed `total_upstream` / `total_downstream` / `total_Tbars`
  (all datasets, including totals-only ones).
* **`<dataset>/roi_mapping.tsv`** — every source ROI / neuropil name seen in the
  data, its `group` (brain/vnc), the VFB region-individual `target` IRI, and how
  many neurons touch it. The group is `unmapped` (and the target blank) for any
  name with no region individual; those are dropped from the OWL. A large
  unmapped primary ROI signals a missing region individual or a broken mapping.
  Written at build time (it needs the source pull), unlike the other reports.

During a build, `report_totals` also prints how each neuron's `total_upstream`
compares with its summed region input. Where the mapped regions partition the
neuron (flywire, banc) the region sum can only be smaller; where they overlap
(hemibrain, male-cns) it can be larger.

**Shape validation** ([`src/sparql/*-violation.sparql`](src/sparql)) locks the
output to exactly the shape we emit:

* region files — a direct `ObjectPropertyAssertion` from a VFB neuron to a VFB
  region individual; every postsynaptic (`RO:0002110`) edge has an `upstream`,
  every presynaptic (`RO:0002113`) edge has a presynaptic count (`Tbars` for
  neuprint, or `downstream` for flywire/banc), and counts are positive; counts
  sit **only** on `upstream`/`Tbars`/`downstream` and **only** on those two
  relations (any `rdf:type`/class-form axiom is rejected); every edge **source**
  is a VFB neuron and every **target** a VFB region individual;
* totals files (`synapse-total-violation.sparql`) — only positive
  `total_upstream` / `total_downstream` / `total_Tbars` assertions, each on a
  VFB neuron.

These are **not** run over the full files (`robot verify` doesn't fit in memory
at connectome scale); every build checks a sample of each file it writes, the
test suite runs them on samples (and asserts the disallowed shapes are
rejected), and the invariants are guaranteed by construction.

> The standard ODK ROBOT report (`<dataset>/obo-report.tsv`, via `make
> obo_reports`) is of limited use here: a side-load edge file intentionally
> carries no labels/definitions for external terms, so the report is almost
> entirely `missing_label`. `obo_reports` filters those rows out.

## Tests & CI

`src/scripts/test_make_region_connectivity.py` covers, with **no live queries**:

* **unit tests** of the pure logic — region-name matching (label prefix/suffix,
  synonyms, ambiguous names, the backslash fix), the BANC neuropil-name mapping,
  edge-count gating, the flywire/banc/catmaid/neuprint loaders (on small local
  files or stubbed responses), totals axioms (including 18-digit root ids),
  reports and header metadata;
* **offline pipeline tests** — network calls and source reads are stubbed with
  *saved* samples in [`src/scripts/test_data/`](src/scripts/test_data), and the
  real transform runs end-to-end through `main()` for neuprint (hemibrain,
  male-cns), flywire and a totals-only dataset; outputs are ROBOT-converted and
  run through the shape checks.

```bash
pip install -r src/scripts/requirements.txt pytest && pytest src/scripts
```

The [`qc.yml`](.github/workflows/qc.yml) GitHub Actions workflow runs this suite
on every push/PR (in the ODK image, no neuprint token needed). It does **not**
rebuild or check the committed artifacts — that needs the live sources and
doesn't fit in a runner's memory.

To refresh the neuprint sample after a neuprint/schema change, regenerate the
files in `test_data/` from a small fixed set of bodyIds.

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
  # totals_only: true                    # totals file only (no roi_dataset needed)
```

The top-level key is the `<dataset>` build arg and the output-file suffix. No
code or Makefile change is needed — `make` picks it up. A neuprint dataset needs
only `site` + `roi_dataset`; a flywire dataset also needs `source: flywire` and
`counts_post_url` / `counts_pre_url`; a banc dataset `source: banc` and
`counts_parquet_url`. A `totals_only` dataset (e.g. `source: catmaid`) needs no
`roi_dataset`; the Makefile builds it with its own rule.

The region individuals must exist in VFB first (loaded from the connectome's
ROIs). Names with no individual are dropped and flagged in `roi_mapping.tsv`.

## Loading

* **Region files** are consumed by
  [`vfb-pipeline-dumps`](https://github.com/VirtualFlyBrain/vfb-pipeline-dumps) as
  edge-only side-loads (kept out of the reasoned/merged graph): its
  `PDB_EXTERNAL_ONTS` list includes `region_connectivity_*.owl`. They are named
  `region_connectivity_*` to distinguish them from the neuron-to-neuron
  `connectome_*` side-loads.
* **Totals files** must go into the **main graph** instead: side-loads are
  converted with `owl2neo4jcsv … only_edges`, which would drop node properties.
  They need listing in
  [`vfb-pipeline-collectdata`](https://github.com/VirtualFlyBrain/vfb-pipeline-collectdata)'s
  `vfb_fullontologies.txt` (e.g.
  `http://virtualflybrain.org/data/VFB/OWL/synapse_totals_*`), like the scRNAseq
  metadata files, and must **not** match `PDB_EXTERNAL_ONTS`.

## Repository layout

| Path | What it is |
|---|---|
| `region_connectivity_<dataset>.owl.gz` | Region-connectivity side-loads (repo root). The uncompressed `.owl` is git-ignored. |
| `synapse_totals_<dataset>.owl.gz` | Whole-neuron synapse totals (node properties), loaded into the main graph. The uncompressed `.owl` is git-ignored. |
| `src/scripts/make_region_connectivity.py` | The generator: source pull → mappings → OWL (region + totals) + ROI report. |
| `src/scripts/cypher_query.py` | Helper for querying VFB's Neo4j (`kb`/`pdb`); also provides the retrying HTTP session. |
| `src/scripts/requirements.txt` | Python dependencies. |
| `src/scripts/test_make_region_connectivity.py` | Unit + offline pipeline tests. |
| `src/scripts/test_data/` | Saved neuprint/kb samples used by the pipeline tests. |
| `src/ontology/dataset_config.yaml` | Per-dataset config — read by the generator and the Makefile. |
| `src/ontology/VFB_region_connectivity.Makefile` | Custom build & report targets (edit this, not `Makefile`). |
| `src/ontology/Makefile`, `run.sh` | ODK-managed wrappers (do not edit). |
| `src/ontology/reports/` | Generated reports — per-dataset in `<dataset>/`, merged term list at the root. |
| `src/sparql/region_connectivity_*.sparql`, `synapse_totals_summary.sparql` | Term-list and summary report queries. |
| `src/sparql/*-violation.sparql` | Shape-validation checks (`synapse-total-violation.sparql` for the totals files). |
| `.github/workflows/qc.yml` | CI: runs the test suite. |

## Acknowledgements

This repository was created using the
[Ontology Development Kit (ODK)](https://github.com/INCATools/ontology-development-kit).

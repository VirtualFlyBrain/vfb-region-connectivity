# SPARQL queries

The queries here run with [`robot`](http://robot.obolibrary.org) over the
generated files (`region_connectivity_<dataset>.owl`,
`synapse_totals_<dataset>.owl`), not over an edited ontology — this repository
has none. They fall into two groups.

## Shape checks (`*-violation.sparql`)

Each returns the rows that break the expected shape of a generated file; any row
fails the check. The build runs them on a sample of every file it writes
(`VFB_region_connectivity.Makefile`), and the test suite runs them on samples
and on deliberately malformed inputs (`src/scripts/test_make_region_connectivity.py`).

| Query | Applies to | Rejects |
|---|---|---|
| `post-edge-missing-upstream` | region files | a `RO:0002110` edge with no `upstream` |
| `pre-edge-missing-count` | region files | a `RO:0002113` edge with neither `Tbars` nor `downstream` |
| `nonpositive-count` | region files | a count <= 0 |
| `unexpected-connectivity-relation` | region files | counts on any relation other than `RO:0002110` / `RO:0002113` (incl. class-form axioms) |
| `unexpected-count-annotation` | region files | any count annotation other than `upstream` / `Tbars` / `downstream` |
| `nonvfb-neuron` | region files | an edge whose source is not a VFB individual |
| `nonvfb-region` | region files | an edge whose target is not a VFB individual |
| `synapse-total-violation` | totals files | a non-positive or unexpected `custom:` total, or one not on a VFB individual |

## Not used

The other queries here (the remaining `*-violation.sparql` files,
`VFB_region_connectivity_terms.sparql`, `basic-report`, `edges`, `terms`,
`synonyms`, `*.ru`, …) are ODK scaffolding for edited ontologies and are not run
by this repository's build.

## Reports

Run by `make region_connectivity_reports`; outputs go to
`src/ontology/reports/<dataset>/` and are committed so versions can be diffed.

| Query | Output |
|---|---|
| `region_connectivity_terms.sparql` | `terms.txt` — neurons and regions in the region file |
| `region_connectivity_summary.sparql` | `summary.tsv` — header provenance, entity/edge counts and count totals of the region file |
| `synapse_totals_summary.sparql` | `totals_summary.tsv` — header provenance, neurons with totals and the summed totals |

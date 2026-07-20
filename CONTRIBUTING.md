# Contributing to VFB Region Connectivity

:+1: Thanks for taking the time to contribute!

This repository is **not a hand-edited ontology** — it is a data pipeline that
pulls per-neuropil synapse counts from [neuprint](https://neuprint.janelia.org)
and generates one side-load OWL file per connectomics dataset. There is no
`-edit` file to edit: the OWL is regenerated from the data. Please read the
[README](README.md) first.

## Reporting a problem

Use the [Issue Tracker](https://github.com/VirtualFlyBrain/vfb-region-connectivity/issues).
The most useful reports for this repo are:

- **A missing or wrong ROI → region mapping** — e.g. an ROI flagged `unmapped` in
  `src/ontology/reports/<dataset>/roi_mapping.tsv` that should map to a VFB region
  individual or an FBbt class (or one mapped to the wrong target).
- **Wrong or unexpected counts** — a neuron's `upstream` / `Tbars` / `downstream`
  in a region that looks incorrect. Note that counts are regenerated live from
  neuprint and drift slightly from older pdb values within the same dataset
  version; that drift is expected, gross discrepancies are not.
- **A new dataset** to add.

See the [issue template](issue_template.md) for the details to include.

## Making a change

The OWL is generated, so fixes are almost always to the pipeline, not the OWL:

- **Mapping / generation logic** → [`src/scripts/make_region_connectivity.py`](src/scripts/make_region_connectivity.py)
- **Build / report targets** → [`src/ontology/VFB_region_connectivity.Makefile`](src/ontology/VFB_region_connectivity.Makefile)
- **Shape-validation checks** → [`src/sparql/`](src/sparql)
- **Adding a dataset** → follow [Adding a dataset](README.md#adding-a-dataset)

Workflow:

1. Fork (external contributors) and create a branch, e.g. `git checkout -b fix-roi-mapping`.
2. Make the change. If you touch the pipeline logic, add/adjust a test in
   `src/scripts/test_make_region_connectivity.py` and run the suite:
   ```bash
   pip install -r src/scripts/requirements.txt pytest && pytest src/scripts
   ```
3. If it affects output, regenerate the affected dataset (`sh run.sh make refresh`)
   and commit the updated `region_connectivity_<dataset>.owl.gz` and the
   `src/ontology/reports/<dataset>/` reports.
4. Open a pull request. CI ([`qc.yml`](.github/workflows/qc.yml)) runs the tests
   and validates the committed artifacts — it must pass.

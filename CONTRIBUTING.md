# Contributing to VFB Region Connectivity

:+1: Thanks for taking the time to contribute!

This repository is **not a hand-edited ontology** — it is a data pipeline that
pulls synapse counts from each connectome's source
([neuprint](https://neuprint.janelia.org), the FlyWire release files on Zenodo,
BANC's public bucket, and VFB's CATMAID servers) and generates, per dataset, a
region-connectivity side-load OWL file and a whole-neuron synapse totals OWL
file. There is no `-edit` file to edit: the OWL is regenerated from the data.
Please read the [README](README.md) first.

## Reporting a problem

Use the [Issue Tracker](https://github.com/VirtualFlyBrain/vfb-region-connectivity/issues).
The most useful reports for this repo are:

- **A missing or wrong ROI → region mapping** — e.g. an ROI flagged `unmapped` in
  `src/ontology/reports/<dataset>/roi_mapping.tsv` that should map to a VFB region
  individual or an FBbt class (or one mapped to the wrong target).
- **Wrong or unexpected counts** — a neuron's `upstream` / `Tbars` / `downstream`
  in a region, or its `total_upstream` / `total_downstream` / `total_Tbars`,
  that looks incorrect. Note that neuprint and CATMAID counts are regenerated
  live and can drift slightly from older pdb values within the same dataset
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
3. If it affects output, regenerate the affected dataset (`sh run.sh make refresh`,
   or delete its `.owl` files and `sh run.sh make`), refresh the reports
   (`sh run.sh make region_connectivity_reports`), and commit the updated
   `region_connectivity_<dataset>.owl.gz`, `synapse_totals_<dataset>.owl.gz` and
   `src/ontology/reports/` files.
4. Open a pull request. CI ([`qc.yml`](.github/workflows/qc.yml)) runs the test
   suite — it must pass. It does not rebuild or check the committed artifacts,
   so check the build's shape-check output and the report diffs yourself.

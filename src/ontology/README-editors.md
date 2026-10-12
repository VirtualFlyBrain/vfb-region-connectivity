These notes are for maintainers of VFB_region_connectivity.

This repository was scaffolded with the
[Ontology Development Kit (ODK)](https://github.com/INCATools/ontology-development-kit),
but it is **not** a hand-edited ontology — there is no `-edit` file and the
standard ODK release targets (`all`, `prepare_release`, `test`) are not used.

Instead it is a data pipeline: it pulls synapse counts from each connectome's
source (neuprint, FlyWire's Zenodo release, BANC's bucket, VFB's CATMAID servers)
and emits, per dataset, a region-connectivity side-load OWL file and a
whole-neuron synapse totals OWL file. Datasets are configured in
[`dataset_config.yaml`](dataset_config.yaml). All of the custom build/report
logic lives in:

* [`VFB_region_connectivity.Makefile`](VFB_region_connectivity.Makefile) — build
  and report targets (`make region_connectivity`, `make region_connectivity_reports`);
* [`../scripts/make_region_connectivity.py`](../scripts/make_region_connectivity.py) — the generator;
* [`../sparql/`](../sparql) — the shape-validation and report queries.

See the [top-level README](../../README.md) for the full description, how to build,
the reports, tests, and how to add a dataset.

These notes are for maintainers of VFB_region_connectivity.

This repository was scaffolded with the
[Ontology Development Kit (ODK)](https://github.com/INCATools/ontology-development-kit),
but it is **not** a hand-edited ontology — there is no `-edit` file and the
standard ODK release targets (`all`, `prepare_release`, `test`) are not used.

Instead it is a data pipeline: it pulls per-neuropil synapse counts from neuprint
and emits one side-load OWL file per connectomics dataset. All of the custom
build/report logic lives in:

* [`VFB_region_connectivity.Makefile`](VFB_region_connectivity.Makefile) — build
  and report targets (`make region_connectivity`, `make region_connectivity_reports`);
* [`../scripts/make_region_connectivity.py`](../scripts/make_region_connectivity.py) — the generator;
* [`../sparql/`](../sparql) — the shape-validation queries.

See the [top-level README](../../README.md) for the full description, how to build,
the reports, tests, and how to add a dataset.

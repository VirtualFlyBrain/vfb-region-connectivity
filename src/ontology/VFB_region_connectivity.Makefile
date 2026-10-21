## Customize Makefile settings for VFB_region_connectivity
##
## This repo does not build a single merged ontology. Instead it produces, per
## connectomics dataset (dataset_config.yaml), up to two OWL files:
##
##     region_connectivity_<dataset>.owl   (written to the repo root)
##     synapse_totals_<dataset>.owl
##
## Each region file holds every neuron's per-neuropil-region synaptic tallies as
## axiom-annotated RO:0002110 / RO:0002113 relations (counts in the n2o
## `custom:` namespace), matching the representation in VFB's pdb. The files are
## side-loaded by vfb-pipeline-dumps (edge-only) and never merged/reasoned, so
## keeping one file per dataset bounds the size of any single artifact.
##
## The same pull also writes synapse_totals_<dataset>.owl: each neuron's
## whole-neuron synapse totals as node annotations. Those must reach pdb as node
## properties, which an edge-only side-load drops, so that file is loaded into
## the main graph instead (vfb-pipeline-collectdata vfb_fullontologies.txt).
## Datasets marked totals_only (CATMAID) produce only the totals file.
##
## Build everything:   sh run.sh make region_connectivity
## Build one dataset:  sh run.sh make ../../region_connectivity_hemibrain.owl
## Totals-only one:    sh run.sh make ../../synapse_totals_catmaid-l1em.owl

# Build the side-load files by default (this repo has no ODK -edit ontology).
.DEFAULT_GOAL := region_connectivity

# Path to a neuprint token (single line), needed by the neuprint datasets.
# See https://neuprint.janelia.org (Account).
NEUPRINT_TOKEN_FILE = ../../neuprint_token.txt

# Datasets to build. Single source of truth is dataset_config.yaml (also read by
# make_region_connectivity.py); DATASETS is its list of top-level keys.
DATASET_CONFIG = dataset_config.yaml
DATASETS := $(shell grep -E '^[A-Za-z0-9._-]+:' $(DATASET_CONFIG) | sed 's/:.*//')
# datasets with `totals_only: true` produce only synapse_totals_<dataset>.owl
TOTALS_ONLY_DATASETS := $(shell awk '/^[A-Za-z0-9._-]+:/{k=$$1; sub(":$$","",k)} /^  totals_only: *true/{print k}' $(DATASET_CONFIG))
REGION_DATASETS := $(filter-out $(TOTALS_ONLY_DATASETS),$(DATASETS))

RELEASEDIR = ../..
REGION_FILES = $(patsubst %, $(RELEASEDIR)/region_connectivity_%.owl, $(REGION_DATASETS))
TOTALS_FILES = $(patsubst %, $(RELEASEDIR)/synapse_totals_%.owl, $(DATASETS))
TOTALS_ONLY_FILES = $(patsubst %, $(RELEASEDIR)/synapse_totals_%.owl, $(TOTALS_ONLY_DATASETS))
# shape checks run on a small QC sample of each generated file (build-time)
RC_QUERIES = $(foreach V,post-edge-missing-upstream pre-edge-missing-count \
	nonpositive-count unexpected-connectivity-relation unexpected-count-annotation \
	nonvfb-neuron nonvfb-region,\
	$(SPARQLDIR)/$(V)-violation.sparql)

######################### Environment #########################

# Python venv for the source pulls. Sentinel rebuilds only when requirements
# change; it is an ORDER-ONLY prereq of the OWL so building reports never
# triggers a (slow) re-pull from the sources.
my-venv/.installed: $(SCRIPTSDIR)/requirements.txt
	python3 -m venv my-venv
	my-venv/bin/pip install -q -r $(SCRIPTSDIR)/requirements.txt
	touch $@

# cypher_query.py is committed; use this to refresh it from the canonical source.
.PHONY: get_query_script
get_query_script:
	wget -O $(SCRIPTSDIR)/cypher_query.py https://raw.githubusercontent.com/VirtualFlyBrain/vfb-scRNAseq-ontology/main/src/scripts/cypher_query.py

######################### Build #########################

.PHONY: region_connectivity
region_connectivity: $(REGION_FILES) $(TOTALS_FILES)

# Force a fresh pull from the sources (data is only rebuilt on demand / when absent).
.PHONY: refresh
refresh:
	rm -f $(REGION_FILES) $(RELEASEDIR)/region_connectivity_*.owl.gz
	rm -f $(TOTALS_FILES) $(RELEASEDIR)/synapse_totals_*.owl.gz
	$(MAKE) region_connectivity

# Build one dataset: pull from its source, emit functional-syntax OWL (one file per
# dataset) plus a small QC sample, convert to RDF/XML + gzip, and run the shape
# checks on the QC sample (the full file is too large for robot verify). The plain
# .owl is the side-load; the .owl.gz is committed. Rebuilt only when the .owl is
# missing (delete it, or `make refresh`, to re-pull). Both files of a dataset
# come from one pull: a multi-target pattern rule runs the recipe once for both.
$(RELEASEDIR)/region_connectivity_%.owl $(RELEASEDIR)/synapse_totals_%.owl: | my-venv/.installed $(TMPDIR)
	my-venv/bin/python3 $(SCRIPTSDIR)/make_region_connectivity.py $* $(NEUPRINT_TOKEN_FILE) $(TMPDIR) $(REPORTDIR)/$*
	$(ROBOT) convert -i $(TMPDIR)/region_connectivity_$*.ofn -f owl -o $(RELEASEDIR)/region_connectivity_$*.owl
	$(ROBOT) convert -i $(TMPDIR)/synapse_totals_$*.ofn -f owl -o $(RELEASEDIR)/synapse_totals_$*.owl
	gzip -9kf $(RELEASEDIR)/region_connectivity_$*.owl $(RELEASEDIR)/synapse_totals_$*.owl
	$(ROBOT) verify -i $(TMPDIR)/qc/region_connectivity_$*.ofn --queries $(RC_QUERIES) -O $(TMPDIR)/qc
	$(ROBOT) verify -i $(TMPDIR)/qc/synapse_totals_$*.ofn --queries $(SPARQLDIR)/synapse-total-violation.sparql -O $(TMPDIR)/qc
	@echo "\nBuilt region_connectivity_$*.owl + synapse_totals_$*.owl (+ .gz), QC passed\n"

# Totals-only datasets (no region file). A static pattern rule, so it takes
# precedence over the two-target pattern rule above for these files.
$(TOTALS_ONLY_FILES): $(RELEASEDIR)/synapse_totals_%.owl: | my-venv/.installed $(TMPDIR)
	my-venv/bin/python3 $(SCRIPTSDIR)/make_region_connectivity.py $* $(NEUPRINT_TOKEN_FILE) $(TMPDIR) $(REPORTDIR)/$*
	$(ROBOT) convert -i $(TMPDIR)/synapse_totals_$*.ofn -f owl -o $@
	gzip -9kf $@
	$(ROBOT) verify -i $(TMPDIR)/qc/synapse_totals_$*.ofn --queries $(SPARQLDIR)/synapse-total-violation.sparql -O $(TMPDIR)/qc
	@echo "\nBuilt $@ (+ .gz), QC passed\n"

######################### Reports #########################
## sh run.sh make region_connectivity_reports
##   -> per-dataset + merged term lists, per-dataset content summaries of the
##   region files, and a per-dataset summary of each synapse_totals file.
##   (shape checks run at build time on a QC sample; the roi_mapping report is
##   written cheaply during the build.)

# Reports load the full OWL, so they are heavy for large datasets; run after a
# build. reports/<dataset>/{terms.txt,summary.tsv}; merged term list at root.
DATASET_TERMS = $(patsubst %, $(REPORTDIR)/%/terms.txt, $(REGION_DATASETS))
DATASET_SUMMARIES = $(patsubst %, $(REPORTDIR)/%/summary.tsv, $(REGION_DATASETS))
TOTALS_SUMMARIES = $(patsubst %, $(REPORTDIR)/%/totals_summary.tsv, $(DATASETS))

$(REPORTDIR)/%/terms.txt: $(RELEASEDIR)/region_connectivity_%.owl | $(REPORTDIR)
	mkdir -p $(REPORTDIR)/$*
	$(ROBOT) query -f csv -i $< --query $(SPARQLDIR)/region_connectivity_terms.sparql $(TMPDIR)/$*-terms.csv
	tail -n +2 $(TMPDIR)/$*-terms.csv | sort -u > $@

# merged term list across all output files
$(REPORTDIR)/region_connectivity-terms.txt: $(DATASET_TERMS) | $(REPORTDIR)
	sort -u $^ > $@

$(REPORTDIR)/%/summary.tsv: $(RELEASEDIR)/region_connectivity_%.owl | $(REPORTDIR)
	mkdir -p $(REPORTDIR)/$*
	$(ROBOT) query -f tsv -i $< --query $(SPARQLDIR)/region_connectivity_summary.sparql $@

# Shape validation is NOT run over the full file: `robot verify` doesn't fit in
# memory at connectome scale (male-cns's ~2 GB OWL OOMs). The same SPARQL checks
# (src/sparql/*-violation.sparql) instead run on a small sample of every generated
# file in src/scripts/test_make_region_connectivity.py, and the invariants are
# guaranteed by construction in build_axioms.

$(REPORTDIR)/%/totals_summary.tsv: $(RELEASEDIR)/synapse_totals_%.owl | $(REPORTDIR)
	mkdir -p $(REPORTDIR)/$*
	$(ROBOT) query -f tsv -i $< --query $(SPARQLDIR)/synapse_totals_summary.sparql $@

.PHONY: region_connectivity_reports
region_connectivity_reports: $(DATASET_TERMS) $(REPORTDIR)/region_connectivity-terms.txt \
	$(DATASET_SUMMARIES) $(TOTALS_SUMMARIES)

# --- optional standard ROBOT report (missing_label rows filtered out: side-load
#     edge files intentionally carry no labels for external terms).
DATASET_OBO_REPORTS = $(patsubst %, $(REPORTDIR)/%/obo-report.tsv, $(REGION_DATASETS))

$(REPORTDIR)/%/obo-report.tsv: $(RELEASEDIR)/region_connectivity_%.owl | $(REPORTDIR)
	mkdir -p $(REPORTDIR)/$*
	$(ROBOT) report -i $< --fail-on none --print 5 -o $(TMPDIR)/$*-obo-report.tsv
	grep -v "	missing_label	" $(TMPDIR)/$*-obo-report.tsv > $@ || cp $(TMPDIR)/$*-obo-report.tsv $@

.PHONY: obo_reports
obo_reports: $(DATASET_OBO_REPORTS)

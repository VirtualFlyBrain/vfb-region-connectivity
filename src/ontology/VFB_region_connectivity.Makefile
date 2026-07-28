## Customize Makefile settings for VFB_region_connectivity
##
## This repo does not build a single merged ontology. Instead it produces one
## side-load OWL file per connectomics dataset:
##
##     region_connectivity_<dataset>.owl   (written to the repo root)
##
## Each file holds every neuron's per-neuropil-region synaptic tallies as
## axiom-annotated RO:0002110 / RO:0002113 relations (counts in the n2o
## `custom:` namespace), matching the representation in VFB's pdb. The files are
## side-loaded by vfb-pipeline-dumps (edge-only) and never merged/reasoned, so
## keeping one file per dataset bounds the size of any single artifact.
##
## Build everything:   sh run.sh make region_connectivity
## Build one dataset:  sh run.sh make ../../region_connectivity_hemibrain.owl

# Build the side-load files by default (this repo has no ODK -edit ontology).
.DEFAULT_GOAL := region_connectivity

# Path to a neuprint token (single line). See https://neuprint.janelia.org (Account).
NEUPRINT_TOKEN_FILE = ../../neuprint_token.txt

# Datasets to build. Single source of truth is dataset_config.yaml (also read by
# make_region_connectivity.py); DATASETS is its list of top-level keys.
DATASET_CONFIG = dataset_config.yaml
DATASETS := $(shell grep -E '^[A-Za-z0-9._-]+:' $(DATASET_CONFIG) | sed 's/:.*//')

RELEASEDIR = ../..
REGION_FILES = $(patsubst %, $(RELEASEDIR)/region_connectivity_%.owl, $(DATASETS))
# shape checks run on a small QC sample of each generated file (build-time)
RC_QUERIES = $(foreach V,post-edge-missing-upstream pre-edge-missing-tbars \
	nonpositive-count unexpected-connectivity-relation unexpected-count-annotation \
	nonvfb-neuron nonvfb-region,\
	$(SPARQLDIR)/$(V)-violation.sparql)

######################### Environment #########################

# Python venv for the neuprint pull. Sentinel rebuilds only when requirements
# change; it is an ORDER-ONLY prereq of the OWL so building reports never
# triggers a (slow) neuprint re-pull.
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
region_connectivity: $(REGION_FILES)

# Force a fresh pull from neuprint (data is only rebuilt on demand / when absent).
.PHONY: refresh
refresh:
	rm -f $(REGION_FILES) $(RELEASEDIR)/region_connectivity_*.owl.gz
	$(MAKE) region_connectivity

# Build one dataset: pull from neuprint, emit functional-syntax OWL (one file per
# dataset) plus a small QC sample, convert to RDF/XML + gzip, and run the shape
# checks on the QC sample (the full file is too large for robot verify). The plain
# .owl is the side-load; the .owl.gz is committed. Rebuilt only when the .owl is
# missing (delete it, or `make refresh`, to re-pull).
$(RELEASEDIR)/region_connectivity_%.owl: | my-venv/.installed $(TMPDIR)
	my-venv/bin/python3 $(SCRIPTSDIR)/make_region_connectivity.py $* $(NEUPRINT_TOKEN_FILE) $(TMPDIR) $(REPORTDIR)/$*
	$(ROBOT) convert -i $(TMPDIR)/region_connectivity_$*.ofn -f owl -o $@
	gzip -kf $@
	$(ROBOT) verify -i $(TMPDIR)/qc/region_connectivity_$*.ofn --queries $(RC_QUERIES) -O $(TMPDIR)/qc
	@echo "\nBuilt $@ (+ $@.gz), QC passed\n"

######################### Reports #########################
## sh run.sh make region_connectivity_reports
##   -> per-dataset + merged term lists, per-dataset content summaries.
##   (shape checks run at build time on a QC sample; the roi_mapping report is
##   written cheaply during the build.)

# Reports load the full OWL, so they are heavy for large datasets; run after a
# build. reports/<dataset>/{terms.txt,summary.tsv}; merged term list at root.
DATASET_TERMS = $(patsubst %, $(REPORTDIR)/%/terms.txt, $(DATASETS))
DATASET_SUMMARIES = $(patsubst %, $(REPORTDIR)/%/summary.tsv, $(DATASETS))

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

.PHONY: region_connectivity_reports
region_connectivity_reports: $(DATASET_TERMS) $(REPORTDIR)/region_connectivity-terms.txt \
	$(DATASET_SUMMARIES)

# --- optional standard ROBOT report (missing_label rows filtered out: side-load
#     edge files intentionally carry no labels for external terms).
DATASET_OBO_REPORTS = $(patsubst %, $(REPORTDIR)/%/obo-report.tsv, $(DATASETS))

$(REPORTDIR)/%/obo-report.tsv: $(RELEASEDIR)/region_connectivity_%.owl | $(REPORTDIR)
	mkdir -p $(REPORTDIR)/$*
	$(ROBOT) report -i $< --fail-on none --print 5 -o $(TMPDIR)/$*-obo-report.tsv
	grep -v "	missing_label	" $(TMPDIR)/$*-obo-report.tsv > $@ || cp $(TMPDIR)/$*-obo-report.tsv $@

.PHONY: obo_reports
obo_reports: $(DATASET_OBO_REPORTS)

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

# Datasets to build (keys must exist in DATASETS in make_region_connectivity.py).
DATASETS = hemibrain

RELEASEDIR = ../..
REGION_FILES = $(patsubst %, $(RELEASEDIR)/region_connectivity_%.owl, $(DATASETS))

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
	rm -f $(REGION_FILES)
	$(MAKE) region_connectivity

# Pull counts from neuprint, emit functional-syntax OWL, normalise to RDF/XML,
# then gzip. The plain .owl (~160 MB for hemibrain) is the side-load artifact;
# the .owl.gz (~3 MB) is what is committed / released (GitHub blocks >100 MB).
# Rebuilt only when missing (delete it, or run `make refresh`, to re-pull).
$(RELEASEDIR)/region_connectivity_%.owl: | my-venv/.installed $(TMPDIR)
	mkdir -p $(REPORTDIR)/$*
	my-venv/bin/python3 $(SCRIPTSDIR)/make_region_connectivity.py $* $(NEUPRINT_TOKEN_FILE) $(TMPDIR)/region_connectivity_$*.ofn $(REPORTDIR)/$*/roi_mapping.tsv
	$(ROBOT) convert -i $(TMPDIR)/region_connectivity_$*.ofn -f owl -o $@
	gzip -kf $@
	@echo "\nBuilt $@ (+ $@.gz)\n"

######################### Reports #########################
## sh run.sh make region_connectivity_reports
##   -> per-dataset + merged term lists, per-dataset content summaries,
##      and shape validation (fails on any violation).

# Dataset-specific reports live in reports/<dataset>/ ; the merged term list sits
# at the reports/ root.

# --- term lists (participating neurons + regions), for version-to-version diffs
DATASET_TERMS = $(patsubst %, $(REPORTDIR)/%/terms.txt, $(DATASETS))

$(REPORTDIR)/%/terms.txt: $(RELEASEDIR)/region_connectivity_%.owl | $(REPORTDIR)
	mkdir -p $(REPORTDIR)/$*
	$(ROBOT) query -f csv -i $< --query $(SPARQLDIR)/region_connectivity_terms.sparql $(TMPDIR)/$*-terms.csv
	tail -n +2 $(TMPDIR)/$*-terms.csv | sort -u > $@

# merged term list across all datasets
$(REPORTDIR)/region_connectivity-terms.txt: $(DATASET_TERMS) | $(REPORTDIR)
	sort -u $^ > $@

# --- content summary (entity/edge counts + count totals), for version diffs
DATASET_SUMMARIES = $(patsubst %, $(REPORTDIR)/%/summary.tsv, $(DATASETS))

$(REPORTDIR)/%/summary.tsv: $(RELEASEDIR)/region_connectivity_%.owl | $(REPORTDIR)
	mkdir -p $(REPORTDIR)/$*
	$(ROBOT) query -f tsv -i $< --query $(SPARQLDIR)/region_connectivity_summary.sparql $@

# --- shape validation: every connectivity edge carries the expected counts, etc.
RC_VALIDATION_CHECKS = post-edge-missing-upstream pre-edge-missing-tbars \
	nonpositive-count unexpected-connectivity-relation nonvfb-neuron
RC_VALIDATION_QUERIES = $(foreach V,$(RC_VALIDATION_CHECKS),$(SPARQLDIR)/$(V)-violation.sparql)
DATASET_VERIFY = $(patsubst %, $(REPORTDIR)/%/verify.txt, $(DATASETS))

$(REPORTDIR)/%/verify.txt: $(RELEASEDIR)/region_connectivity_%.owl | $(REPORTDIR)
	mkdir -p $(REPORTDIR)/$*
	$(ROBOT) verify -i $< --queries $(RC_VALIDATION_QUERIES) -O $(REPORTDIR)/$*
	echo "PASS" > $@

.PHONY: region_connectivity_reports
region_connectivity_reports: $(DATASET_TERMS) $(REPORTDIR)/region_connectivity-terms.txt \
	$(DATASET_SUMMARIES) $(DATASET_VERIFY)

# --- optional standard ROBOT report (missing_label rows filtered out: side-load
#     edge files intentionally carry no labels for external terms).
DATASET_OBO_REPORTS = $(patsubst %, $(REPORTDIR)/%/obo-report.tsv, $(DATASETS))

$(REPORTDIR)/%/obo-report.tsv: $(RELEASEDIR)/region_connectivity_%.owl | $(REPORTDIR)
	mkdir -p $(REPORTDIR)/$*
	$(ROBOT) report -i $< --fail-on none --print 5 -o $(TMPDIR)/$*-obo-report.tsv
	grep -v "	missing_label	" $(TMPDIR)/$*-obo-report.tsv > $@ || cp $(TMPDIR)/$*-obo-report.tsv $@

.PHONY: obo_reports
obo_reports: $(DATASET_OBO_REPORTS)

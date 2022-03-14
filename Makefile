# This Makefile builds things for all analysis code at once. It supports the
# following targets:
#   - all: Run `make` in every analysis folder
#   - clean: Run `make clean` in every analysis folder
#   - Pruners: Run `make Pruners` in every analysis folder
#   - <Analysis>: Builds everything related to the given analysis
# The Makefile will always try to compile the `pruner` executable in Bonzais/.
#
# WARNING: This Makefile won't fail if some analysis code is broken (ie doesn't
#          compile). It's not a bug, it's a feature.

# List of the analysis to compile
ANALYSIS = DYJets HZZ2l2nu TagAndProbe #WJets

# Error out on old make
MIN_VERSION := 4.0
OK := $(filter $(MIN_VERSION),$(firstword $(sort $(MAKE_VERSION) $(MIN_VERSION))))
ifeq ("x$(OK)", "x")
$(error 'The version of make you are using ($(MAKE_VERSION)) is too old. Run cmsenv')
endif

# This rule builds everything
ALL = $(addsuffix .all, $(ANALYSIS))
.PHONY: all
all: Bonzais.all $(ALL)
%.all: Bonzais.all
	$(MAKE) -C $(basename $@) || true

# This rule cleans everything
CLEAN = $(addsuffix .clean, $(ANALYSIS))
.PHONY: clean
clean: $(CLEAN) Bonzais.clean
%.clean:
	$(MAKE) -C $(basename $@) clean || true

# This rule builds all pruners
PRUNERS = $(addsuffix .pruners, $(ANALYSIS))
.PHONY: Pruners
Pruners: $(PRUNERS)
%.pruners: Bonzais.all
	$(MAKE) -C $(basename $@) Pruners || true

# Analysis rules (eg `make DYJets`)
$(foreach ANA,$(ANALYSIS),$(eval .PHONY: $(ANA)))
$(foreach ANA,$(ANALYSIS),$(eval $(ANA): $(ANA).all))

Welcome to Shears, the Simple and Handy Event Analysis ROOT-based Suite
=======================================================================

Shears is the analysis framework developed within the CMS SMP-VJ group for Run-2 data analyses. It is an evolution of the W+jets and Z+jets analysis framework of 8 TeV analysis [1] originally developed by ULB (Tomislav Seva and Alexandre Leonard) and including contributions from METU (Bugra Bilin), Northeastern University (Apichart Hortiangthan) and CEA/Saclay (Philippe Gras).

More information can be found in the dedicated [Shears Twiki](https://twiki.cern.ch/twiki/bin/view/CMS/ShearsAnalysisFramework): 

Repository organisation
-----------------------

**Bonzais** Tools to make small ROOT trees ("skims") from genuine trees (so called Boababs)

**DYJets** Analysis code. Currently for W+jet and Z+jet analysis.

**HZZ2l2nu** Analysis code from the HZZ2l2nu group.

**TagAndProbe** Analysis code for tag and probe (also from the HZZ2l2nu group).

**WJets** W+jet analysis code (can run on both Baobab and Bonzai formats, currently set for Bonzai).

**ntuple_production** Tools to produce the Boabab ROOT ntuple.

**Baobabs** The code that produces the ntuple (Boabab) from CMS EDM dataset (MINIAOD) and used by the tools in ntuple\_production

*Documentation can be found in the respective directories (README.md file displayed when browsing the directory with gitlab as for the one you are reading).*

Installation recipe
------------------

### If you don't want to produce baobabs

First, download the code (change the URL if don't use the main shears):

```
git clone ssh://git@gitlab.cern.ch:7999/shears/shears.git
```

Once you have the code, you can checkout another branch as usual. For example, in order to run on 2017 data:

```
cd shears
git checkout Run2017
cd ..
```

You can then compile (see [below](#compiling)).

### If you want to baobabs

First, download the code (change the URL if don't use the main shears):

```
git clone ssh://git@gitlab.cern.ch:7999/shears/shears.git
```

Once you have the code, you can checkout another branch as usual. For example, in order to run on 2017 data:

```
cd shears
git checkout Run2017
cd ..
```

There is a script in the `shears/Baobabs` directory that takes care of installing the required CMSSW packages. Run it:

```
shears/Baobabs/install.sh
```

The script will:

* Checkout the required CMSSW version in the current directory
* Install all required modules
* Move shears into the CMSSW tree (in `src/`) and replace it by a soft link
* Compile everything

If you don't need it, you can remove the soft link:

```
rm shears
```

Environment
-----------

Shears needs a working (and decently recent) CMSSW environment. Find a CMSSW installation and set up the environment using:

```
cd your_cmssw_directory
cmsenv
cd -
```

Shears provides some scripts under the `shears/ntuple_production` directory. You can add them to your `PATH` by using:

```
cd shears
export PATH="$PATH:$(readlink -e ntuple_production)"
```

Compiling
---------

[Setup the environment](#environment) and use `make` to compile the code.

In addition to analysis' own `Makefile`s, the framework provides a `Makefile` at
the root of the source tree. It supports the following targets:

~~~{.sh}
make            # Will build everything for all analysis. Probably not what you want
make all        # Same as above
make clean      # Will clean everything (calls `make clean` in analysis folders)
make Pruners    # Will build all pruners (by running `make Pruners` in every
                # analysis folder)
make <Analysis> # Will build everything related to the given analysis (ie the pruner
                # executable, `make Pruners` and `make all` in the analysis folder)
~~~

Other targets should be considered *internal* and not relied on (althrough they
may work).

Old code
--------

Old code is pruned off the repo from time to time, but is still available in the
history. The table below lists the last commit at which such features were
available.

| Commit   | Comment                         |
|----------|---------------------------------|
| a0a6677b | Code for reading 8 TeV baobabs  |

References
----------

[1] https://github.com/iihe-cms-sw/TreeAnalysis

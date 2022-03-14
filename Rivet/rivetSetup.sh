#!/bin/bash
if [ -d $CMSSW_BASE/src/GeneratorInterface/RivetInterface/test/rivetSetup.sh ]; then
    source $CMSSW_BASE/src/GeneratorInterface/RivetInterface/test/rivetSetup.sh
else
    source $CMSSW_RELEASE_BASE/src/GeneratorInterface/RivetInterface/test/rivetSetup.sh
fi


for GROUP in Top FSQ SMP HIG SUS HIN Exotica BPH B2G
do
  export RIVET_REF_PATH=$RIVET_REF_PATH:$CMSSW_BASE/src/Rivet/${GROUP}/data
  export RIVET_INFO_PATH=$RIVET_INFO_PATH:$CMSSW_BASE/src/Rivet/${GROUP}/data
  export RIVET_PLOT_PATH=$RIVET_PLOT_PATH:$CMSSW_BASE/src/Rivet/${GROUP}/data
done

which yodamerge &> /dev/null || GETYODA=1
if [ $GETYODA==1 ]; then
  eval `scram tool info yoda | grep YODA_BASE`
  export PATH=$PATH:$YODA_BASE/bin
fi

# cmsRivet scripts
export PATH=$PATH:$CMSSW_BASE/src/Rivet/scripts

#!/bin/bash
set -x
set -e

mybasename="`basename "$0"`"
myfullpath="`readlink -f "$0"`"
mydir="`dirname "$myfullpath"`"

function set_exit_code(){
# At the end of the script modify the FJR
    exitCode="90001"
    errorType="$1"
    exitMessage="$2"
    if [ -e FrameworkJobReport.xml ]
    then
	cat <<EOF >FrameworkJobReport.xml.tmp
<FrameworkJobReport>
<FrameworkError ExitStatus="$exitCode" Type="$errorType" >
$exitMessage
</FrameworkError>
EOF
	tail -n+2 FrameworkJobReport.xml >> FrameworkJobReport.xml.tmp
	mv FrameworkJobReport.xml.tmp FrameworkJobReport.xml
    else
	cat <<EOF > FrameworkJobReport.xml
<FrameworkJobReport>
<FrameworkError ExitStatus="$exitCode" Type="$errorType" >
$exitMessage
</FrameworkError>
</FrameworkJobReport>
EOF
    fi
}

die(){
    echo "$@" 1>&2
    set_exit_code 1 scriptError "$*"
    exit 1
}


date;
t1=`date +%s`

#echo Arguments:
#echo "$@"
#echo "----------------------------------------------------------------------"
#echo
#
#echo Environment:
#env
#echo "----------------------------------------------------------------------"
#echo
#
#echo 'Working directory content (xdev, max-depth 3):'
#find . -xdev -maxdepth 3
#echo "----------------------------------------------------------------------"
#echo

#produces FramworkJobReport.xml
#It's a dummy run and we don't need to loop
#on any event.
cat > myPSet.py <<EOF
import FWCore.ParameterSet.Config as cms
import pickle
process = pickle.load(open('PSet.pkl', 'rb'))
process.maxEvents.input = 0
EOF
cmsRun -j FrameworkJobReport.xml myPSet.py

NJob="$1"

if [ "$NJob" = 0 ]; then
    echo "NJob=0! Forced to 1"
    NJob=1
fi

shift

while [ $# -gt 0 ]; do
    echo "$1" | grep -q '^[[:alpha:]_][[:alnum:]_]*=\([^[:space:]]*\|"[^"]*"\)$'
    [ $? = 0 ] && eval "$1"
    shift
done

echo "cfg=$cfg"
echo "maxEvents=$maxEvents"
[ -n "$cfg" ] || die "Parameter cfg was not found!"

unset maxEventsOpt
[ -n "$maxEvents" ] && maxEventsOpt="maxEvents=$maxEvents"

[ -n "$NJob" ] || die "Missing job ID"

echo "Job id: $NJob"

#note: when using --dryrun option of crab submit, the job is run twice in the same directory, we therefore
#need to look for libRooUnfold.so both in local directory and RooUnfold one, where it is moved to by this
#script.
[ -f libRooUnfold.so -o -f RooUnfold/libRooUnfold.so ] || die "You need to add libRooUnfold.so in your task input file list."

mkdir RooUnfold
mv libRooUnfold.so RooUnfold/
mv RooUnfoldDict_rdict.pcm RooUnfold/


tar xzf EfficiencyTables.tgz
tar xzf rcdata.2016.v3.tgz

#%lep% keyword in the is used to provide to configurations, on for DMu and one for DE
echo "$cfg" | grep -q lepSel  && lepSels="DMu" || lepSels="dummy"

nRuns=20
#if [ $NJob -gt $nRuns ]; then
#    iRun=$((NJob-nRuns))
#    lepSel=DE
#else
iRun=$NJob
lepSel=DMu
#fi

export VJETS_CONFIG="`echo "$cfg" | sed "s/lepSel/${lepSel}/"`"
echo "Running with configuraion file $VJETS_CONFIG..."
    
case "$iRun" in
    1)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=1 ;;
    2)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=2 ;;
    3)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=3 ;;
    4)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=4 ;;
    5)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=5 ;;
    6)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=6 ;;
    7)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=7 ;;
    8)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=8 ;;
    9)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=9 ;;
    10)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=10 ;;
    11)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=11 ;;
    12)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=12 ;;
    13)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=13 ;;
    14)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=14 ;;
    15)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=15 ;;
    16)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=16 ;;
    17)  ./runZJets_newformat $maxEventsOpt doWhat=DYJETS whichSyst=0 nJobs=17 jobNum=17 ;;

esac

tar czf HistoFiles.tgz HistoFiles*

echo "List of files:"
ls


date 
t2=`date +%s`
echo "Duration:  $(((t2*10-t1*10+300)/600)) mn"


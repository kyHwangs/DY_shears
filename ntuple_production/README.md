Scripts to produce Baobab and Bonzai Ntuples
===========================================

Author: Ph. Gras. CEA/Saclay

Setting up working environment of ntuple producton
=================================================

1. If you don't have yet a copy of the Shears framework, download it from the [gitlab repository](https://gitlab.cern.ch/shears/shears):
```
git clone ssh://git@gitlab.cern.ch:7999/shears/shears.git
```
2. Set crab 3 environment following instructions from <https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3CheatSheet#Environment_setup>. On lxplus for bash shell user:
```
source /cvmfs/cms.cern.ch/crab3/crab_light.sh
```
3. Set up Tuple environment following instructions from <https://github.com/UGent/Tupel/blob/Tupel_MiniAOD/Tupel/README.md>.
4. Add the `shears/ntuple_production` to your command search path PATH. Assuming you use the bash shell:
```
cd shears/ntuple_production
PATH=$PATH:`pwd`
cd -
```

Important: if you don't use the "light" version of the crab 3 environmemt, the CMSSW environment must be set after the crab 3 environment. Otherwise you will inherit from a Python version too old for the Baobab production tools.

Production of Baobab ntuples
============================

You should first follow the instruction of the "Setting up environment..." section. The command to steer the Baobab ntuple production is `grow_boabab`. The data sets to process, the json file in case of real data and the EOS destination directory for the produce ntuple need to be listed in a file in the format defined in [1], we suggest to call datasets.txt. Follows an example of datasets.txt file content:

```
# output directory: /store/group/phys_smp/AnalysisFramework/Test/Ntuple
# json file: /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_*-*_13TeV_PromptReco_Collisions15_25ns_JSON.txt

/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD
```

The output directory path must end with the command

`grow_baobabs --set-catalog-path datasets.txt`

To submit the jobs to produce the Boabab ntuple of the datasets listed in the catalog file, runs:

`grow_baobabs --new-jobs`  # if there is an error like: urllib2.HTTPError: HTTP Error 403: Forbidden---- need to do "voms-proxy-init -voms cms"

To check running job and generate luminosity section summary json file, use the --check option:

`grow_baobabs --check`

The production is tracked in the baobab\_prod.sqlite3` database file. For a production, it is important to run grow\_baobabs always in the same disk directory as this database file is picked up form the working directory.

To list the content of the production database file, used the --list option:

`grow_baobabs --list`

A [Crab 3](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCrab) directory is created for each GRID task (processing of a dataset, called also <it>job batch</i>). It is name crab\_ followed by the dataset name and the job batch id. Crab 3 commands can be used directly. Use the Crab option -d to specify the task directory.

For central production operation, you should also read the instructions specific to the central production which can be in, https://twiki.cern.ch/twiki/bin/view/CMS/SmpVjBaobabProduction.

More options for grow\_baobabs command
--------------------------------------

Some usual options: 

--no-submit, to be used with --new-jobs to only creates the crab configuration files and register the crab tasks in the grow\_baobabs database. The job can then be submitted in a second step with the --submit command or using directly the CRAB3 crab command;

--task, to limit command to a given task: see in the online help to which commands it applies;

--resubmit, to resubmit jobs;

--no-crab-report-for-mc, to be used with the --check command to disable the crab report retrieval for MC dataset, which can take ages;

More options are available, run `grow_baobabs --help` to get the online help.

The simple\_grow\_baobabs command
-------------------------------

The command simple\_grow\_baobabs is a simplified version of grow\_baobab, which can be used for private production. It does not keep track of the production as grow\_baobabs does and does have a built-in support for the integrated luminosity calucation. For a central production, grow\_baobabs should be used.

**Beware**: the destination directory, where the Baobab ntuples will be written to, is specified in the catalog file. If you copy a central production dataset.txt file for a private production, don't forget to edit the file and modify the *output directory* line.


The mgrow\_boababs command
-------------------------

This command will run grow\_baobabs with the provided arguments in each directory defined in the prodrc configuration file (dirs variable). Online help can be obtained by running `warden --help`

The warden command
------------------

Warden is a command to monitor a Shear ntuple production. Online help can be obtained by running `warden --help`

The list of monitored directory is defined in the prodrc file located in the same directory than warden.

Note 1: when the warden command is run, it creates two files, warden_summary and warden_last_log in the directory it is executing from.

Note 2: if warden is running in background (loop mode), you can stop independently under which login it is running, with the warden --stop command. The command will work on two conditions: the warden command is located in the same directory that when it was started and you have write access to this directory.

The getlumi command
-------------------

This command will retrieve the integrated luminosity corresponding to a list of luminosity section provided in a json file. It is wrapper to CMS standard commands and is used by the grow\_boababs command.  It will be updated to follow latest recommendation for integrated luminosity calculation. Online help can be obtained by running `getlumi --help`.

The new-era command
-------------------

This command is used to set up work space on AFS and EOS for SMP Shears ntuple production. Online help can be obtained by running `new-era --help`.


Production of Bonzai ntuples
============================

Bonzai ntuple can be produced using the [Pruner](https://gitlab.cern.ch/shears/shears/tree/master/Bonzais/Pruner) utility. The grow_bonzai tool can be used for massing production using the GRID infractuction. The grow_bonzai tool is less advanced that the grow_baobab ones. Both tools will be eventually merged.


Usage:
-----

* You should first follow the instruction of the "Setting up environment..." section.
* Create a clean working directory and enter into this directory
```mkdir bonzai-prod
cd bonzai-prod
```
* Create a task list file like the [grow_bonzai_task_list_example.txt](grow_bonzai_task_list_example.txt). You should specify in the file the EOS location to write the ntuple to.
* Run:

```grow_bonzais --task-list your_task_list.txt```

* Use the standard crab command to check the GRID task status or use the script `get_crab_status`. To use this script you should remove from your directory old crab_XXXX subdirectories and keep only the ones of your current submissions. This script will ask you for each process if the submission is terminated. If you answer yes, the crab task will be recorded in the jobs_ok file and the task won't be check at the next call to `get_crab_status`. Once you are done remove the file jobs_ok, otherwise if your submit again tasks with the same names from the same directory, the `get_crab_status` will ignore them.

* Once the jobs are successfully completed, you can produce the catalog files using the command:

```grow_bonzais --make-catalogs your_task_list.txt```

If you want to produce the catalog of completed jobs before all jobs are done, you can comment the uncompleted ones in the task list file. For instance add at the beginning of the line the comment %running% .

Online help can be obtained by executing `grow_bonzai --help`. As for grow_baobab, it is important to set up CMSSW environment after CRAB3.

Common Utilities
================

The cmslog-inputfiles command
----------------------------

This is a utility to retrieve the list of input files used by a crab job from its log file. Online help can be obtained by running `cmslog-inputfile --help`.

The dumpPSetpkl command
-----------------------

This is a utility to dump the CMSSW configuration store by Crab in  input/PSet.pkl. Online help can be obtained by running `dumpPSetpkl --help`.

Other files
----------

send_shears_emails: this script is used by the warden command to send announcement of new ntuple to the cms-shears-ntuple-announcement egroup.

prodrc: configuration file for warden and mgrow\_boababs

[1] <https://twiki.cern.ch/twiki/bin/view/CMS/SmpVjBaobabProduction>

Submission of V+jet analysis on the grid
========================================

This directory contains a list of tool to submit the production of V+jet analysis histograms on the grid using crab.

Files:
------

### Provided files:

<dl>
<dt><code>crabConfig.py</code>:</dt>    <dd>crab configuration file</dd>
<dt><code>job_crab.sh</code>:</dt>      <dd>shell script which is run on the grid nodes</dd>
<dt><code>donothing_cfg.py</code>:</dt> <dd>dummy CMSSW configuration required by the job</dd>
</dl>

### Additionnal file to be provided:

* The `runZJets_newformat` configuration file(s)
* A tarball with efficiency tables (see instructions below)

Instructions
------------

To submit on the grid, you to fist set a CMSSW environment and the crab environment. The CMSSW version/architecture combination must be suppored by the Grid. If crab complains it is not the case, try another combination (*).

1. If you use another environment than you use to, recompile all the code to make sure it is compiled with the compiler version from the new CMSSW environment. Run `make clean` in the `DYJets` directory, then run `make` in `RooUnfold` and finally run `make -j 16` in the `DYJets` directory (`-j 16` argument is optional).
2. Copy your configuration file in the current directory. The default job runs on both muon and electron channels. You should provide a file for each channel. Include DE and DMu in the names of the respective files, e.&nbsp;g. `vjet_DE.cfg` and `vjet_DMu.cfg`.
3. Produce a tarball with the efficiency tables: in the `DYJets` directory, run `tar czf EfficiencyTables.tgz EfficiencyTable`
4. Edit the crabConfig.py file to set the parameters `config.JobType.inputFiles` and `config.Data.outLFNDirBase`. You can also change the parameter `config.Site.storageSite` if you want to write the output on another site than CERN.
5. Submit your tasks by running the command `crab submit`
6. Please refer to the crab documentation for information on how to monitor the jobs.

Output
------

Each job will create a tarball that will be written in the output directory you set in the crab configuration file (`config.Data.outLFNDirBase`). The tarball contains the `histoDir` directory.

Customising the jobs
--------------------

The default script, `job_crab.sh`, runs the Z+jet analysis for all the standard samples and all the systematic uncertainties for both electron and muon channels. To run with different options, you will need to edit the job_crab.sh script. The list of job command lines can be found in the `case` block. It is important to:

* change the variable nRuns value in the script accordingly to the number of entries in the case block;
* change the number of units (`config.Data.totalUnit`) in the `crabConfig.py` file to matches the number of jobs: the number of entries in the case block times the number of channels.

Processing a sample in pieces
-----------------------------

For large sample, like the signal MC samples, it can be useful to split the sample in pieces processed by different jobs. This can be achieved by using the `nJobs` and `jobNum` options of `runZJets_newfomat` program. An example of crab configuration and script can be found in this directory (`DYJets/crab.d`): `crabConfig.py`, `job_crabDY.sh` and for the analysis configuration file `vjets_silver_lepSel_crabDY.cfg`. The number of samples pieces is defined by the nJobs variable of the script, you can find around line 118. If its value is changed the number of total number of jobs should be updated in the crab configuration file. Note that this number is displayed at the beginning of the log of the jobs to ease the check of the crab configuration. You can check the number a posteriori or run a local test with the command `./job_crabDY.sh 1 cfg=vjets_silver_lepSel_crabDY.cfg maxEvents=1` to let the script to compute it for you.

Note
----

If you use the `--dryrun` option of `crab submit` to test the jobs, beware it will run on the full set of data events. The reason is that the `job_crab.sh` script does not use the maximum number of events set by crab in the CMSSW configuration file. You add temporally in your crab configuration the option `maxEvents=10` in the list of arguments defined by `config.JobType.scriptArgs` to limit the run to 10 events.

job_crab.sh script and add a `maxEvents` option to limit the number of events.

(*) I don't know how to get the list of supported combinations. If you are aware of how to get it, I'm interested to learn it. -- Philippe.

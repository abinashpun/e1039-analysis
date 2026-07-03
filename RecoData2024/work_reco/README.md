# e1039-analysis/RecoData2024/work_reco

The spill-level data files are stored in `/pnfs/e1039/persistent/users/kenichi/dst`.
These data are used for reconstruction, as selected by `DIR_DST` in `gridsub.sh`.


## Test of Reconstruction Using a Few Runs

The procedure below will analyze the 1st good spill (not run) to draw basic dimuon variables.

```
cd /path/to/your/e1039-analysis/RecoData2024/work_reco
source ../setup.sh
./gridsub.sh -j 1-1
root -b 'AnaEventTree.C("scratch/reco")'
display result/h1_m.png &
```

Here are what happens in each step:
- `gridsub.sh` sets up a subdirectory (per run) to analyze DST files of the given run(s).
    - The option `-j 1-1` selects the 1st spill (listed in `list_run_spill.txt`).
    - Processing one spill is directed by `gridrun.sh` and `Fun4All.C`.
    - A SubsysReco module, `AnaDimuon`, is registered by default.
      Thus a simple analysis is carried out after the reconstruction.
    - The process runs on local computer by default (i.e. when `-g` is not given).
    - Dimuon parameters extracted by `AnaDimuon` are stored in `scratch/reco/run_*/out/output.root` as TTree.
- `AnaEventTree.C` reads the TTree files to draw histograms.
    - Histograms are saved under `result/`.
    - Histograms should be empty since only one spill was processed so far.


## Reconstruction of All Runs/Spills

You can process all runs/spills by adjusting the `-j` option of `gridsub.sh`.
The `-g` option let you use the Grid computing.
You first submit jobs for the first 10 spills;
```
kinit
./gridsub.sh -g -j 1-10
```

You might execute `jobsub_q_mine` to monitor the job status, which is an alias of `jobsub_q --group spinquest --user=$USER` defined in `setup.sh`.
Once you confirm that the jobs are running fine, you submit more jobs for the remaining runs;
```
./gridsub.sh -g -m 150 -j 11-
```
Here `-m 150` sets the max number of total (idle+running+held) jobs at a time.
This command takes a long time (like a night) to finish, so you had better execute it on `screen` or in background.

The output files are stored under `data/reco/` in case Grid is used.
You can analyze them once (almost) all jobs finish;
```
root -b AnaEventTree.C
display result/h1_m.png &
```


## Reconstruction Condition

The reconstruction condition is programmed in `Fun4All.C`.
Thus you need look into the macro to understand the default condition.
You might modify the macro as you need.


## Re-Analysis of Reconstructed Data

The analysis defined in `AnaDimuon` is done together with the reconstruction.
When you need reanalyzed the reconstructed data, you can use `work_ana/`.

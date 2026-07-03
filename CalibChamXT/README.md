# e1039-analysis/CalibChamXT

A program to calibrate the chamber X-T curve, using cosmic-ray events at present.
The file structure, the configuration and the usage are similar to the following programs;

* https://cdcvs.fnal.gov/redmine/projects/calib-cham-rt
* https://github.com/E1039-Collaboration/e1039-analysis/tree/master/AnaRealDst
* https://github.com/E1039-Collaboration/e1039-analysis/tree/master/RecoDev
* /seaquest/users/apun/abi_project/data_manage/e1039-data-mgt_test


## Download

You are recommended to use `spinquestgpvm01`.
You first clone the GitHub repository as usual;

```
ssh spinquestgpvm01.fnal.gov
cd /path/to/your/working/directory
git clone git@github.com:E1039-Collaboration/e1039-analysis.git
cd e1039-analysis/CalibChamXT
```

The top directory of this program is `e1039-analysis/CalibChamXT/`.
All files/directories explained below are relative to the top directory, unless specified.


## Preparation

The runs to be analyzed are listed in `work/list_run.txt`.
The list contains only run 3492 at present.
You could modify the list when necessary.

You copy the DST files of all runs from `data2` to `pnfs` as follows.
It usually takes a few minutes per run;

```
cd work
./copy_dst.sh
```

The source code is located in `src/`.
The following commands set up the shell environment and compile the code;

```
source setup.sh
cmake-this
make-this
```


## Track Reconstruction

The program can execute the track reconstruction, not yet the calibration of the X-T curve.

The following commands carry out a test, using 10,000 events in the 1st DST file.
It should take 10-20 seconds.

```
cd work
./run_reco.sh -d 1 -e 10000 GO
```

The main output file is `scratch/main/run_003492_spill_000000000_spin/out/DSTreco.root`.
The ROOT file (`DSTreco.root`) contains all tracklets.
You had better check the number of tree entries is non-zero.

You then execute the following commands to process all events/runs.
The `-g` option is to run jobs on the grid.

```
kinit
./run_reco.sh -g GO
```

You use `jobsub_q --group spinquest --user=$USER` (or `./check_job.sh`) to wait until all grid jobs finish.

```
~kenichi/bin/hadddraw -o hist_raw data/main/run_*/out/output_*.root
```


## Tracklet Analysis

Then the following commands analyze the tracklets reconstructed by `run_reco.sh`.
Plots will appear in `calib/1/`.

```
cd work
./make_rt.sh
```

This script executes an Fun4All macro (`Fun4AllAna.C`), which makes use of a SubsysReco modle (`src/SRMakeRTCurve`) to analyze `DSTreco.root`.

Note that `make_rt.sh` has two commented-out lines that find and analyze the evaluation file (`eval.root`) using a ROOT macro (`AnaRTCurve.C`).
This part is kept for backward compatibility for now, but probably be deleted soon.


## To-Do List

1. Clean up old sources/macros/scripts that use the evaluation file.
1. Implement an iterative extraction of X-T curves.


## Reference

Station ID:  https://github.com/E1039-Collaboration/e1039-wiki/wiki/event-reconstruction

Trigger bits:
* FPGA1: Any ST2 and Any ST4
* FPGA2: Top ST2 and Top ST4
* FPGA3: Top ST2 and Bottom ST4
* FPGA4: Bottom ST2 and Top ST4
* FPGA5: ((H1X||H1Y)&&(H2X||H2Y)) || ((H2X||H2Y)&&(H4X||H4Y2))
* NIM1: (H1X||H1Y)&&(H2X||H2Y)&&(H3X||H4Y1)&&(H4X||H4Y2)
* NIM2: (H1X||H1Y)&&(H2X||H2Y)
* NIM3: Random
* NIM4: (H2X||H2Y)&&(H4X||H4Y2)

## Author

Kenichi Naknao <knakano0524@gmail.com>

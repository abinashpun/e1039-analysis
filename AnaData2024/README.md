# e1039-analysis/AnaData2024

A sample program to analyze the commissioning data.
It analyzes the good runs and spills listed in `work/list_run_spill.txt`, which were defined in [DocDB 10854-v2](https://seaquest-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=10854&version=2).
It reads the spill-level DST files stored in `/pnfs/e1039/persistent/users/kenichi/dst`.


## Test of Analysis Using a Few Runs

The procedure below will analyze the first three good runs to extract chamber hit distributions.

```
# Compile the source code.
cd /path/to/your/e1039-analysis/AnaData2024
source setup.sh
cmake-this
make-this

# Extract chamber hits from DST to TTree.
cd work
./gridsub.sh -j 1-3

# Draw hit distributions.
root -b 'AnaEventTree.C("scratch/default")'
display result/h1_ele_30.png &
```

Here are what happens in each step:
- `source setup.sh` does
    - Set up the E1039 software environment, and 
    - Define shell variables and functions to build and use the source code.
    - You need execute this command when opening a new text terminal.
- `cmake-this` and `make-this` build the source code under `src/`.
    You can execute these commands at any directory.
- `gridsub.sh` sets up a subdirectory per run to analyze DST files of the runs.
    - The option `-j 1-3` selects the 1st-3rd runs (found in `list_run_spill.txt`).
    - Analyzing one run is directed by `gridrun.sh` and `Fun4All.C`.
    - A SubsysReco module, `AnaChamHit`, is registered by default, to extract chamber hits.
    - The analysis runs on local computer by default (i.e. when `-g` is not given).
    - Extracted chamber hits are stored in `scratch/default/run_*/out/output.root` as TTree.
- `AnaEventTree.C` reads the TTree files to draw histograms.
    - Histograms are saved under `result/`.


## Analysis Using All Runs

You can analyze other runs by adjusting the `-j` option of `gridsub.sh`.
The `-g` option let you use the Grid computing.
You first submit jobs for the first 10 runs;
```
cd work
kinit
./gridsub.sh -g -j 1-10
```

You might execute `jobsub_q_mine` to monitor the job status, which is an alias of `jobsub_q --group spinquest --user=$USER` defined in `setup.sh`.
Once you confirm that the jobs are running fine, you submit more jobs for the remaining runs;
```
./gridsub.sh -g -j 11-
```

The output files are stored under `data/` in case Grid is used.
You can analyze them once (almost) all jobs finish;
```
root -b AnaEventTree.C
display result/h1_ele_30.png &
```


## Modification of Analysis

The contents of analysis are programmed in `src/AnaChamHit.(h|cc)`.
Thus you need modify these files to do the analysis you need.
You need execute `make-this` after modification.

Or you might create a new SubsysReco module (like `src/AnaHodoHit.(h|cc)`) and register it in `Fun4All.C`.
You need execute `cmake-this` after creating a new module and `make-this` after modifying it.


## Contact

Kenichi Naknao <knakano0524@gmail.com>

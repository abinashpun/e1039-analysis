# e1039-analysis/RecalRealDst

A compact program to recalibrate the real data (in the DST format).


## Basic Usage

The procedure below will
- Read the DST files of runs 6155 and 6156 under `/pnfs/e1039/persistent/users/kenichi/dst`,
- Calibrate the element position and the drift distance of hits, and
- Save the recalibrated DST files under `/pnfs/e1039/scratch/users/$USER/recal`.

```
cd /path/to/your/e1039-analysis/AnaRealDst
source setup.sh
cmake-this
make-this
cd work
./make-input-list.sh
./do-recal.sh
```

You should modify `make-input-list.sh` to add more runs to be processed.
You then execute `./make-input-list.sh` and `./do-recal.sh` again, where `do-recal.sh` skip the runs (spills) whose output file already exists.


## Advanced Usage

- The directory of input DST files is defined with `DIR_IN` in `make-input-list.sh`.
  You might change it.
- The directory of output DST files is defined with `DIR_OUT_BASE` in `do-recal.sh`.
  You might change it.
- The calibrations applied are selected (registered) in `Fun4All.C`.
- You can create your own calibrator under `src/` and register it in `Fun4All.C`.
  `CalibDriftDistMod` is an example.


## Contact

Kenichi Naknao <knakano0524@gmail.com>

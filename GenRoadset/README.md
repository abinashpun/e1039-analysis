# e1039-analysis/GenRoadset

A program to generate a roadset (i.e. FPGA trigger logic) for E1039,
using signal event and background (BG) event.
This package was made from `knakano0524/e1039-trigger-tune`.

By default, the Drell-Yan simulation is used to generate signal events, and
the E906-NIM3 data are used as BG events.
The full-BG simulation will become usable to generate BG events.

## Source Files

All source files are placed in `src/`.
You use the following commands to compile them.

```
source setup.sh
cmake-this
make-this
```


## Analysis Procedure

1. Generation of signal events
1. Generation of BG events
1. Generation of new roadset
1. Evaluation of new roadset

The 1st step is carried out in `macro_gen_signal/`.
You can find details in [`macro_gen_signal/README.md`](macro_gen_signal/README.md).

The 2nd step is carried out in `macro_e906_nim3/`.
You can find details in [`macro_e906_nim3/README.md`](macro_e906_nim3/README.md).
But you usually need not carry out this step, since a full set of extracted events is already available.
The existing set is used in the next step by default.

The 3rd and 4th steps are carried out in `macro_ana/`.
You can find details in [`macro_ana/README.md`](macro_ana/README.md).


## To-Do List

* Take and use RF ID & Turn ID for NIM3 events.


## Contact

Kenichi Nakano <knakano0524@gmail.com>

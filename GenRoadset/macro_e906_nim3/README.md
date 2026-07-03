# Macro to extract E906 NIM3 events

This directory is to extract NIM3 events from the E906 digit files.
It outputs a small ROOT file named `bg_tree.root`, which will be used later in `macro_ana`.

You usually need not carry out this step by yourself.


## Setting of Event Extraction

The setting of the event extraction is configured in `Fun4All.C` as usual.
The SubsysReco module that creates `bg_tree.root` is `src/SubsysRecoBG.(h|cc)`.
You can modify them as you need.


## Procedure

The digit files had better be prestaged before use.

```
../../HitEmbedding/macro_gen_emb_e906/prestage-input-file.sh
```

The generation process can/should be done on local computer, since one run takes only a few minutes.
It might be better to divide the whole process into multiple blocks like this;

``
./gridsub.sh -n 1-100
./gridsub.sh -n 101-1000
./gridsub.sh -n 1001-2000
./gridsub.sh -n 2001-4000
``

The run-by-run output files should be merged;

```
./merge-output-file.sh
```

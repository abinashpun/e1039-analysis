# Macro to generate BG events from the full-BG simulation

This directory is to generate BG events from the full-BG simulation.
It outputs a small ROOT file named `bg_tree.root`, which will be used later in `macro_ana`.


## Setting of Event Generation

The setting of the event generation is configured in `Fun4All.C` as usual.
The SubsysReco module that creates `bg_tree.root` is `src/SubsysRecoBG.(h|cc)`.
You can modify them as you need.

The configuration of `SQPileupGen` was taken from `abinashpun/e1039-abi-projects/PileupDev` on 2021-12-21.
`e906_rf00.root` was copied from it also.


## Test Production

You first execute the following command to execute the macro on local for test,
which runs one job that generates 10 in-acceptance events;

```
./gridsub.sh
```


## Mass Production

You then execute the following commands to submit grid jobs (`-g`),
which runs 100 jobs (`-j 100`) that generate events as many as possible (`-e 0`);

```
kinit
./gridsub.sh -g -j 100 -e 0
```

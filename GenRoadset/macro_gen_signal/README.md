# Macro to generate signal events

This directory is to generate signal events to be used for the roadset tuning.
It outputs a small ROOT file named `signal_tree.root`, which will be used later in `macro_ana`.


## Setting of Event Generation

The setting of the event generation is configured in `Fun4All.C` as usual.
The SubsysReco module that creates `signal_tree.root` is `src/SubsysRecoSignal.(h|cc)`.
You can modify them as you need.


## Test Production

You first execute the following command to execute the macro on local for test,
which runs one job that generates 10 in-acceptance events;

```
./gridsub.sh
```


## Mass Production

You then execute the following commands to submit grid jobs (`-g`),
which runs 100 jobs (`-j 100`) that generate 1000 events/job (`-e 1000`);

```
kinit
./gridsub.sh -g -j 100 -e 1000
```

You might submit more jobs when the statistics are not enough;

```
./gridsub.sh -g -j 101-200 -e 1000
```

We usually submit 800 jobs to have the sufficient statistics, as of May 2024.

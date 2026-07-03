# Macro for Event Generation

This directory is to generate signal events.
The usage is similar to `e1039-analysis/SimChainDev`.
You generate sets of J/psi and psi' events at least,
and optionally another set of Drell-Yan events to draw the mass distributions.
The outputs are `DST.root`, which will be used later in `macro_ana/`.


## Setting of Event Generation

The setting of the event generation is configured in `Fun4All.C` as usual.
But note that the number of events/job that you set is of _accepted_ events, not _all_generated_ events.
The followings are recommended numbers of events/job and jobs:

* `SQPrimaryParticleGen`
    - Drell-Yan: 500 events/job = 0.5-3 hours, 200 jobs
    - J/psi    : 500 events/job = 0.5-3 hours, 200 jobs
    - psi'     : 500 events/job = 0.5-3 hours, 200 jobs
* PYTHIA8
    - Drell-Yan: to be written.
    - J/psi    : to be written.
    - psi'     : to be written.


## Test Production

You first execute the following command to execute the macro on local for test,
which runs one job that generates 10 in-acceptance events;

```
./gridsub.sh
```


## Mass Production

You then execute the following commands to submit grid jobs (`-g`),
which runs 100 jobs (`-j 100`) that generate 500 events/job (`-e 500`);

```
kinit
./gridsub.sh -g -j 200 -e 500
```

You might submit more jobs when the statistics are not enough;

```
./gridsub.sh -g -j 201-400 -e 500
```

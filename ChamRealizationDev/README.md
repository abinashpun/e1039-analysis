# e1039-analysis/ChamRealizationDev

A compact program to develop the chamber realization, which adds inefficiencies and position resolutions to chamber and prop tube hits.


## Preparation

The version of e1039-core is selected in `setup.sh`.
The pre-built offline version is used by default.
If you need modify e1039-core, you build e1039-core and select it in `setup.sh`.


## Source Build

The source code is located in `src/`.
The following commands set up the shell environment and compile the code;

```
source setup.sh
cmake-this
make-this
```

`AnaTrack` is a SubsysReco module to analyze reconstructed tracks.
The track reconstruction and the analysis are done before and after the chamber realization in order to evalaute the effect of the chamber realization.


## Event Generation

```
cd macro
root -b 'Fun4AllStep1.C(10000)'
```

This step generates a set of single-track events using `PHG4SimpleEventGenerator`.
The vertex and the momentum of tracks are fixed by default for simplicity.
You can change them as you need, by modifying `Fun4AllStep1.C`.

This step also includes the reconstruction and the analysis of tracks.
Results of the analysis will appear under `AnaTrack_raw/`.
The chamber realization is carried out in the next step.


## Chamber Realization

```
root -b Fun4AllStep2.C
```

This step reads the DST file (`DST.root`) created in the previous step.
It applies the chamber realization and then reconstruct and analyze tracks.
Results of the analysis will appear under `AnaTrack_real/`.

You compare `AnaTrack_raw/` and `AnaTrack_real/` to evaluate the effect of the chamber realization.


## Contact

Kenichi Naknao <knakano0524@gmail.com>

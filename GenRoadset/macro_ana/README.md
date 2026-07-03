# Macro to generate and evaluate roadset

This directory is to first generate a new roadset, using the signal and BG events.
It then evaluates the new roadset, namely re-analyzes the signal and BG events to compute the signal acceptance and the BG rate.

You need execute `source setup.sh` beforehand.


## Lists of Input Files

Signal files and BG files should be listed in `list_signal.txt` and `list_bg.txt`.
You can make the lists by the following command;

```
./make-input-list.sh
```

When you use the full-BG simulation, you have to change `DIR_IN_BG` in the script beforehand.

If you have modified some setting in `macro_gen_signal` or `macro_e906_nim3`, you might have to modify `make-input-list.sh` accordingly.


## QIE Inhibit Threshold

You set the QIE inhibit threshold when generating and evaluating a new roadset.
It was 1200 on average in the last E906 run, although it varied with time.
Thus you should try 1200 and also different values (like 1000 and 1400) to check how the result changes.


## Generation of Roadset

Typical procedure:

```
root -b GenerateRoad.C
```

Related files:
* `macro_ana/GenerateRoad.C`
* `src/AnaGmc.(h|cc)`
* `src/AnaNim3.(h|cc)`
* `src/GenRoadset.(h|cc)`


## Evaluation of Roadset

Typical procedure:

```
root -b ApplyRoadset2Signal.C
root -b 'ApplyRoadset2BG.C(0)'
root -b 'ApplyRoadset2BG.C(45000)'
root -b 'ApplyRoadset2BG.C(63000)'
root -b 'ApplyRoadset2BG.C(80500)'
root -b DrawRS.C
```

Related files:
* `work/ApplyRoadset2Signal.C`
* `work/ApplyRoasetd2BG.C`
* `work/DrawRS.C`
* `src/ReAnaSignal.(h|cc)`
* `src/ReAnaBG.(h|cc)`
* `src/DrawRoadset.(h|cc)`


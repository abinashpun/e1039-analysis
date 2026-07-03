# work_dist

This directory is to extract the X-T curve using the time-distribution method, without the track reconstruction.
The outcome is typically used as the initial X-T curve in the calibration.

## Dataset

By default the macros in this directory uses run 3824 for D0.

## Usage

```
root -b 'Fun4All.C(3824, "/data2/e1039/dst/run_003824_spin.root")'
root -b ExtractXTCurve.C
```

The 1st step takes 17 minutes when all events are analyzed.

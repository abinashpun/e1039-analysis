# e1039-analysis/DstTreeIO

A test program to read and write the "DST" file.
It does not use the Fun4All framework to 
* Build the data structure nor
* Loop over events.

This scheme might be useful in making random (not sequential) access to events.
But I don't yet recommend to use it for serious analyses 
because I don't find whether/how we apply calibrations/selections, which are implemented and shared as `SubsysReco` modules.


## Reading DST

```
source setup.sh
root -b ReadDstTree.C
```

It reads the DST file of run 4621 and prints out basic variables.
It does not create any output file.


## Writing DST

```
source setup.sh
root -b WriteDstTree.C
```

It creates `DSTx.root` and fill the SQRun, SQEvent, SQHitVector and SQTriggerHitVector nodes with dummy values.
The created file is able to be processed with `AnaRealDst` for example.


## Contact

Kenichi Naknao <knakano0524@gmail.com>

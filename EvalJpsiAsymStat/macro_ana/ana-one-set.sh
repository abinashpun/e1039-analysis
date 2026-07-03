#!/bin/bash
PROC=$1 # jpsi, psip or dy
rm -rf   $PROC
mkdir -p $PROC
cd       $PROC
../make-dst-list.sh ../../macro_gen/data/$PROC
root -b -q ../Fun4All.C
root -b -q "../AnaTree.C(\"$PROC\")"

#!/bin/bash
## A temporary script to execute 'MakeRT.C'.
## You need wait for all jobs (submitted by 'run_reco.sh') finished
## before executing this script.
DIR_BASE=$(dirname $(readlink -f $BASH_SOURCE))

ITER=1
LIST_RUN="list_run.txt"
VERSION="main"
DIR_DATA=/pnfs/e1039/scratch/users/$USER/CalibChamXT/$VERSION

LIST_DATA=scratch/list_input_data.txt
mkdir -p scratch

RUN1=0 # 1st run number
while read RUN ; do
    test $RUN1 -eq 0 && RUN1=$RUN
    RUN6=$(printf '%06d' $RUN)
    find $DIR_DATA/run_${RUN6}_* -name DSTreco.root
    #find $DIR_DATA/run_${RUN6}_* -name eval.root
done <$LIST_RUN >$LIST_DATA

# Use the 1st run number, assuming the condition stays the same.
root.exe -b -q "Fun4AllAna.C($RUN1, $ITER, \"$LIST_DATA\")"
#root.exe -b -q "AnaRTCurve.C($RUN1, $ITER, \"$LIST_DATA\")"

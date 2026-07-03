#!/bin/bash

DIR_LOG=log
FN_LIST=list_input.txt
DIR_OUT_BASE=/pnfs/e1039/scratch/users/$USER/recal
#DIR_OUT_BASE=/data4/e1039_data/users/$USER/recal
N_PROC=2

function ProcessOneFile {
    local -r     II=$1
    local -r    RUN=$2
    local -r  SPILL=$3
    local -r  DST_IN=$4
    local -r   RUN6=$(printf '%06d' $RUN)
    local -r SPILL9=$(printf '%09d' $SPILL)
    local -r DIR_OUT=$DIR_OUT_BASE/run_$RUN6
    local -r  FN_OUT=$DIR_OUT/run_${RUN6}_spill_${SPILL9}_spin.root # DST
    #local -r  FN_OUT=$DIR_OUT/run_${RUN6}_spill_${SPILL9}_sraw.root # SRawEvent
    local -r  FN_LOG=$DIR_LOG/run_${RUN6}_spill_${SPILL9}_log.txt
    printf "%6d / %d: run %d  spill %d" $II $NN $RUN $SPILL # $DST $FN_OUT"
    if [ -e $FN_OUT ] ; then
	echo ' --- Skipped as already exists.'
	return
    fi
    echo 
    mkdir -p $DIR_OUT
    root.exe -b -q "Fun4All.C($RUN, $SPILL, \"$DST_IN\", \"$FN_OUT\")" &>$FN_LOG
}

function WaitByNumProcs {
    local -r N_MAX=$1
    while true ; do
        #echo "----------------"
        N_JOBS=$(jobs -r | wc -l)
        echo "N of procs = $N_JOBS"
        test $N_JOBS -le $N_MAX && break
        sleep 4
    done
}

##
## Main
##
mkdir -p $DIR_LOG
NN=$(cat $FN_LIST | wc -l)
II=0
while read RUN SPILL DST ; do
    (( II++ ))
    WaitByNumProcs $N_PROC
    ProcessOneFile $II $RUN $SPILL $DST &
    sleep 2
done <$FN_LIST
WaitByNumProcs 0

#!/bin/bash
DIR_IN=/pnfs/e1039/persistent/users/kenichi/dst # /data2/e1039/dst
FN_OUT=list_input.txt

function ParseOneRun {
    local -r DIR_IN=$1
    local -r RUN=$2
    local -r RUN6=$(printf '%06d' $RUN)
    test -d $DIR_IN/run_$RUN6 || return
    for FNAME in $DIR_IN/run_$RUN6/run_${RUN6}_spill_*_spin.root ; do
	BNAME=$(basename $FNAME)
	SPILL9=${BNAME:17:9}
	SPILL=$(( 10#$SPILL9 ))
    	echo "$RUN $SPILL $FNAME"
    done
}

## In case you process a list of runs
for RUN in 6155 6156 ; do
    echo "Run $RUN" >/dev/stderr
    ParseOneRun $DIR_IN $RUN
done >$FN_OUT

## In case you process a range of runs
#for (( RUN = 6156 ; RUN >= 5846 ; RUN-- )) ; do
#    echo "Run $RUN" >/dev/stderr
#    ParseOneRun $DIR_IN $RUN
#done >$FN_OUT


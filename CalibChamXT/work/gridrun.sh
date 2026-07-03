#!/bin/bash

RUN=$1
ITER=$2
FN_DST=$3
N_EVT=$4

if [ -z "$CONDOR_DIR_INPUT" -o -z "$CONDOR_DIR_OUTPUT" ] ; then
    echo "!ERROR!  CONDOR_DIR_INPUT/OUTPUT is undefined.  Abort."
    exit 1
fi
echo "INPUT  = $CONDOR_DIR_INPUT"
echo "OUTPUT = $CONDOR_DIR_OUTPUT"
echo "RUN    = $RUN"
echo "ITER   = $ITER"
echo "FN_DST = $FN_DST"
echo "N_EVT  = $N_EVT"
echo "HOST   = $HOSTNAME"
echo "PWD    = $PWD"

tar xzf $CONDOR_DIR_INPUT/input.tar.gz
source setup.sh
echo "E1039_CORE = $E1039_CORE"

if [ ${FN_DST:0:1} != '/' ] ; then # in case of relative path
    FN_DST=$CONDOR_DIR_INPUT/$FN_DST
fi

touch timestamp.txt
time root.exe -b -q "Fun4AllReco.C($RUN, $ITER, \"$FN_DST\", $N_EVT)"
RET=$?
if [ $RET -ne 0 ] ; then
    echo "Error in Fun4AllReco.C: $RET"
    exit $RET
fi

find . -mindepth 1 -maxdepth 1 -newer timestamp.txt -exec mv {} $CONDOR_DIR_OUTPUT \;

echo "gridrun.sh finished!"

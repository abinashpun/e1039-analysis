#!/bin/bash

N_EVT=$1
FN_SIG=DST.root
FN_EMB=embedding_data.root

if [ -z "$CONDOR_DIR_INPUT" -o -z "$CONDOR_DIR_OUTPUT" ] ; then
    echo "!ERROR!  CONDOR_DIR_INPUT/OUTPUT is undefined.  Abort."
    exit 1
fi
echo "INPUT  = $CONDOR_DIR_INPUT"
echo "OUTPUT = $CONDOR_DIR_OUTPUT"
echo "HOST   = $HOSTNAME"
echo "PWD    = $PWD"

tar xzf $CONDOR_DIR_INPUT/input.tar.gz
source setup.sh

touch timestamp.txt # All files created after this will be moved to CONDOR_DIR_OUTPUT

time root -b -q "Fun4Sim.C(\"$CONDOR_DIR_INPUT/$FN_SIG\", \"$CONDOR_DIR_INPUT/$FN_EMB\", $N_EVT)"
RET=$?
if [ $RET -ne 0 ] ; then
    echo "Error in Fun4Sim.C: $RET"
    exit $RET
fi

echo "$RET" >status.txt

find . -mindepth 1 -maxdepth 1 -newer timestamp.txt -exec mv {} $CONDOR_DIR_OUTPUT \;

echo "gridrun.sh finished!"

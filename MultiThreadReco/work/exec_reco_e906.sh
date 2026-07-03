#!/bin/bash
DIR_SCRIPT=$(dirname $(readlink -f $0))
DIR_DAT=/data2/users/kenichi/2024-03-27-split-data/output

# Run 28693, Spill 1415125-1415195
# Run 28694, Spill 1415209-1415238
RUN=${1:-'28693'}
SPILL=${2:-'1415125'}
N_EVT=${3:-'0'}

source $DIR_SCRIPT/../setup.sh

RUN6=$(printf '%06d' $RUN)
SPILL9=$(printf '%09d' $SPILL)

FN_IN=$DIR_DAT/digit_run_${RUN6}_spill_${SPILL9}.root

DIR_OUT=$DIR_SCRIPT/out
FN_OUT=$DIR_OUT/srec_run_${RUN6}_spill_${SPILL9}.root

DIR_LOG=$DIR_SCRIPT/log
FN_LOG=$DIR_LOG/log_run_${RUN6}_spill_${SPILL9}.txt

umask 0002
mkdir -p $DIR_OUT $DIR_LOG
{
    TIMEFORMAT='TIME: real %R user %U sys %S'
    time root.exe -b -q "$DIR_SCRIPT/RecoE906Data.C($RUN, $SPILL, \"$FN_IN\", \"$FN_OUT\", $N_EVT)"
    echo "RET: $?"
} |& tee $FN_LOG
#exit $RET

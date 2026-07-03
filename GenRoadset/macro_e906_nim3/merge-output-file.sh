#!/bin/bash
## Script to merge the output files, bg_data.root.
DIR_MACRO=$(dirname $(readlink -f $BASH_SOURCE))

RUN_LIST=/data2/production/list/R009/R009_fy2017.list
DIR_IN=$DIR_MACRO/scratch
DIR_OUT=/pnfs/e1039/persistent/users/$USER/GenRoadsetCommon
FN_OUT=bg_data.root

if [ -e $DIR_OUT ] ; then
    echo "DIR_OUT exists."
    if [ "X$1" != 'Xgo' ] ; then
	echo "Give 'go' to clean it up and proceed.  Abort."
	exit
    fi
    echo "  Clean up..."
    rm -rf $DIR_OUT
fi

##
## List up input files
##
source /e906/app/software/osg/software/e1039/this-e1039.sh
N_RUN=0
LIST_ROOT=
while read DB_SERVER DB_SCHEMA RUN ROADSET ; do
    RUN6=$(printf '%06d' $RUN)
    DIR_RUN=$DIR_IN/${RUN6:0:2}/${RUN6:2:2}/${RUN6:4:2}/out
    FN_RUN=$DIR_RUN/bg_data.root
    test -e $FN_RUN || continue
    
    (( N_RUN++ ))
    printf "Run %6d: %4d\n" $RUN $N_RUN
    LIST_ROOT+=" $FN_RUN"
done <$RUN_LIST

##
## Run hadd
##
echo "Write to $DIR_OUT/$FN_OUT..."
mkdir -p $DIR_OUT
hadd -k $DIR_OUT/$FN_OUT $LIST_ROOT &>$DIR_OUT/log_hadd.txt
    
{
    echo "DIR_MACRO = $DIR_MACRO"
    echo "RUN_LIST  = $RUN_LIST"
    echo "DIR_IN    = $DIR_IN"
    echo "DIR_OUT   = $DIR_OUT"
    echo "FN_OUT    = $FN_OUT"
    echo "N_RUN     = $N_RUN"
    date
} | tee $DIR_OUT/info.txt

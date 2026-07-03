#!/bin/bash
## Script to copy DST files from '/data2' to '/pnfs'.

LIST_RUN="list_run.txt"
DIR_SRC=/data2/e1039/dst
DIR_DEST=/pnfs/e1039/scratch/users/$USER/CalibChamXT/dst
echo "LIST_RUN = $LIST_RUN"
echo "DIR_SRC  = $DIR_SRC"
echo "DIR_DEST = $DIR_DEST"

while read RUN ; do
    echo "Run $RUN"
    RUN6=$(printf '%06d' $RUN)
    if [ -e $DIR_DEST/run_$RUN6 ] ; then
	echo "  Already exist in DIR_DEST.  Skipped."
	continue
    fi
    echo "  Copying..."
    mkdir -p $DIR_DEST
    cp -a $DIR_SRC/run_$RUN6 $DIR_DEST
    echo "  ...done."
done <$LIST_RUN

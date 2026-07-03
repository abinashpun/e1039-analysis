#!/bin/bash
## A temporary script to check the status of jobs submitted.
DIR_BASE=$(dirname $(readlink -f $BASH_SOURCE))

DIR_DATA=/pnfs/e1039/scratch/users/$USER/CalibChamXT/main
DIR_TMP=scratch
LOOP=no

OPTIND=1
while getopts ":d:l" OPT ; do
    case $OPT in
        d ) DIR_DATA=$OPTARG ;;
        l ) LOOP=yes ;;
    esac
done
shift $((OPTIND - 1))


if [ ! -e $DIR_DATA/input.tar.gz ] ; then
    echo "Cannot find '$DIR_DATA/input.tar.gz'.  Not yet submit jobs?  Abort."
    exit
fi

mkdir -p    $DIR_TMP
LIST_JOB_ID=$DIR_TMP/list_job_id.txt
echo "LIST_JOB_ID = $LIST_JOB_ID"
if [ ! -e $LIST_JOB_ID -o $LIST_JOB_ID -ot $DIR_DATA/input.tar.gz ] ; then
    echo "Updating the job-ID list..."
    grep '^Use job id ' $DIR_DATA/*/log_jobsub_submit.txt | awk '{print $4}' >$LIST_JOB_ID
    echo "...done."
fi

N_JOB_ALL=$(cat $LIST_JOB_ID | wc -l)
echo "N_JOB_ALL  = $N_JOB_ALL"

while true ; do
    N_JOB_IDLE=0
    N_JOB_RUN=0
    N_JOB_HELD=0
    while read ID OWNER DATE TIME RUN ST PRI SIZE CMD ; do
	case $ST in
	    I ) (( N_JOB_IDLE++ )) ;;
	    R ) (( N_JOB_RUN++ )) ;;
	    H ) (( N_JOB_HELD++ )) ;;
	esac
    done < <(jobsub_q --group spinquest --user=$USER | grep -f $LIST_JOB_ID)
    echo "N_JOB_IDLE = $N_JOB_IDLE"
    echo "N_JOB_RUN  = $N_JOB_RUN"
    echo "N_JOB_HELD = $N_JOB_HELD"

    test $LOOP = no && break
    test $N_JOB_IDLE -eq 0 -a $N_JOB_RUN -eq 0 && break
    echo "---- Wait ------------------------------------------------------------"
    sleep 300
done

exit 0

#!/bin/bash
DIR_MACRO=$(dirname $(readlink -f $BASH_SOURCE))

#DIR_BKG="/pnfs/e1039/persistent/users/apun/bkg_study/e1039pythiaGen_17Nov21"
#DIR_BKG="/pnfs/e1039/persistent/users/apun/bkg_study/e1039pythiaGen_Jan22/003"
DIR_BKG="/pnfs/e1039/persistent/users/apun/bkg_study/fullbg_candidates/01"

#KMAG_POL=+1
#JOB_NAME=main_01_nim3potv2_H1Xgap # Normal KMag polarity

#KMAG_POL=-1
#JOB_NAME=main_01_nim3pot_reverseMag_H1Xgap
# Reverse KMag polarity #main_01_run6nim3_pcoeff025_thes70k #main_jan22_003

#KMAG_POL=0
#JOB_NAME=main_01_nim3pot_zero_H1Xgap

KMAG_POL=+1
KMAG_SC=0.5
JOB_NAME=main_01_nim3pot_normal050

DO_OVERWRITE=no
USE_GRID=no
JOB_B=1
JOB_E=1
N_EVT=10
OPTIND=1
while getopts ":n:ogj:e:" OPT ; do
    case $OPT in
	n ) JOB_NAME=$OPTARG ;;
	o ) DO_OVERWRITE=yes ;;
        g ) USE_GRID=yes ;;
        j ) JOB_E=$OPTARG ;;
        e ) N_EVT=$OPTARG ;;
    esac
done
shift $((OPTIND - 1))

if [ "${JOB_E%-*}" != "$JOB_E" ] ; then # Contain '-'
    JOB_B=${JOB_E%-*} # Before '-'
    JOB_E=${JOB_E#*-} # After '-'
fi

echo "JOB_NAME     = $JOB_NAME"
echo "DO_OVERWRITE = $DO_OVERWRITE"
echo "USE_GRID     = $USE_GRID"
echo "JOB_B...E    = $JOB_B...$JOB_E"
echo "N_EVT        = $N_EVT"
if [ $USE_GRID == yes ]; then
    DIR_DATA=/pnfs/e1039/scratch/users/$USER/GenRoadset/data_full_bg
    DIR_WORK=$DIR_DATA/$JOB_NAME
    ln -nfs $DIR_DATA data # for convenience
else
    DIR_WORK=$DIR_MACRO/scratch/$JOB_NAME
fi

cd $DIR_MACRO
mkdir -p $DIR_WORK
rm -f    $DIR_WORK/input.tar.gz
tar czf  $DIR_WORK/input.tar.gz  *.C ../inst e906_rf00*.root nim3pot_run6.root

for (( JOB_I = $JOB_B; JOB_I <= $JOB_E; JOB_I++ )) ; do
    DIR_WORK_JOB=$DIR_WORK/$(printf "%04d" $JOB_I)
    if [ -e $DIR_WORK_JOB ] ; then
	echo -n "  DIR_WORK_JOB already exists."
	if [ $DO_OVERWRITE = yes ] ; then
	    echo "  Clean up."
	    rm -rf $DIR_WORK_JOB
	else
	    echo "  Skip."
	    continue
	fi
    fi

    #FN_BKG="${JOB_I}_bkge1039_pythia_17Nov21_100M.root"
    #FN_BKG="${JOB_I}_bkge1039_pythia_Jan22_100M.root"
   
    #runid=`expr $JOB_I + 1000`
    #FN_BKG="0${runid}_bkg_100M_v1.root"  

    FN_BKG=$(printf '%05i_bkg_100M_v1.root' $JOB_I)
 
    mkdir -p $DIR_WORK_JOB/out
    cp -p $DIR_MACRO/gridrun.sh $DIR_WORK_JOB
    
    if [ $USE_GRID == yes ]; then
	CMD="/exp/seaquest/app/software/script/jobsub_submit_spinquest.sh"
	CMD+=" --expected-lifetime='medium'" # medium=8h, short=3h, long=23h
	#CMD+=" --memory=1GB"
	#CMD+=" --lines '+FERMIHTC_AutoRelease=True'"
	#CMD+=" --lines '+FERMIHTC_GraceLifetime=7200' --lines '+FERMIHTC_GraceMemory=2048'" #2 hours of grace lifetime and 2 GB of grace memory if hold
	CMD+=" -L $DIR_WORK_JOB/log_gridrun.txt"
	CMD+=" -f $DIR_WORK/input.tar.gz"
	CMD+=" -f $DIR_BKG/$FN_BKG"
	CMD+=" -d OUTPUT $DIR_WORK_JOB/out"
	CMD+=" file://$DIR_WORK_JOB/gridrun.sh $FN_BKG $N_EVT $KMAG_POL $KMAG_SC"
	unbuffer $CMD |& tee $DIR_WORK_JOB/log_jobsub_submit.txt
	RET_SUB=${PIPESTATUS[0]}
	test $RET_SUB -ne 0 && exit $RET_SUB
    else
	export  CONDOR_DIR_INPUT=$DIR_WORK_JOB/in
	export CONDOR_DIR_OUTPUT=$DIR_WORK_JOB/out
	mkdir -p $DIR_WORK_JOB/in
	cp -p $DIR_WORK/input.tar.gz $DIR_WORK_JOB/in
	ln -nfs $DIR_BKG/$FN_BKG $DIR_WORK_JOB/in
	cd $DIR_WORK_JOB
	$DIR_WORK_JOB/gridrun.sh $FN_BKG $N_EVT $KMAG_POL $KMAG_SC |& tee $DIR_WORK_JOB/log_gridrun.txt
    fi
done

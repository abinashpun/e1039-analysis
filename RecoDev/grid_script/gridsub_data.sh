#!/bin/bash

dir_macros=$(dirname $(readlink -f $BASH_SOURCE))
LIFE_TIME=long # short (3h), medium (8h) or long (23h)

jobname=$1
do_sub=$2

run_list=$3
nevents=$4

echo "njobs=$njobs"
echo "nevents=$nevents"

 
if [ $do_sub == 1 ]; then
    echo "Grid mode."
    if ! which jobsub_submit &>/dev/null ; then
	echo "Command 'jobsub_submit' not found."
	echo "Forget 'source /e906/app/software/script/setup-jobsub-spinquest.sh'?"
	exit
    fi
    work=/pnfs/e906/persistent/users/$USER/RecoDataDev/$jobname
    ln -sf /pnfs/e906/persistent/users/$USER/RecoDataDev data
else
    echo "Local mode."
    work=$dir_macros/scratch/$jobname
fi

mkdir -p $work
chmod -R 01755 $work

cd $dir_macros
tar -czvf $work/input.tar.gz *.C *root
cd -


for run_num in $(cat $run_list) 
do  
 
#location of the data in your persistent area
data_dir="/pnfs/e906/persistent/users/$USER/DstRun"
data_file=$(printf 'run_%06d_spin.root' $run_num)

#use this format for the splitted runs
#data_dir="location/to/your/splitted/long_run"
#data_file=$(printf 'run_001283_spill_%09d_spin.root' $run_num) #for 1283 long run


  mkdir -p $work/$run_num/log
  mkdir -p $work/$run_num/out
  chmod -R 01755 $work/$run_num

  rsync -av $dir_macros/gridrun_data.sh $work/$run_num/gridrun_data.sh

  if [ $do_sub == 1 ]; then
    cmd="jobsub_submit"
    cmd="$cmd -g --OS=SL7 --use_gftp --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE -e IFDHC_VERSION --expected-lifetime='$LIFE_TIME'"
    cmd="$cmd --mail_never"
    cmd="$cmd -L $work/$run_num/log/log.txt"
    cmd="$cmd -f $work/input.tar.gz"
    cmd="$cmd -d OUTPUT $work/$run_num/out"
    cmd="$cmd --append_condor_requirements='(TARGET.GLIDEIN_Site isnt \"UCSD\")'"
    cmd="$cmd -f $data_dir/$data_file"
    cmd="$cmd file://`which $work/$run_num/gridrun_data.sh` $nevents $run_num"
  
    echo "$cmd"
    $cmd
  else
    mkdir -p $work/$run_num/input
    rsync -av $work/input.tar.gz $work/$run_num/input
    cd $work/$run_num/
    $work/$run_num/gridrun_data.sh $nevents $run_num | tee $work/$run_num/log/log.txt
    cd -
  fi
done 2>&1 | tee log_gridsub.txt

## When your job fails due to bad grid nodes,
## you can use the following option to exclude those nodes;
##   cmd="$cmd --append_condor_requirements='(TARGET.GLIDEIN_Site isnt \"UCSD\")'"
## Valid site names are listed here;
## https://cdcvs.fnal.gov/redmine/projects/fife/wiki/Information_about_job_submission_to_OSG_sites
## According to the Fermilab Service Desk, the "--blacklist" option has a known defect.

#!/bin/bash

dir_macros=$(dirname $(readlink -f $BASH_SOURCE))

jobname=$1
do_sub=$2
run_list=$3
nevents=$4
dst_mode=${5:-'split'} # 'split' or 'single'

echo "njobs=$njobs"
echo "nevents=$nevents"

 
if [ $do_sub == 1 ]; then
    echo "Grid mode."
    work=/pnfs/e1039/scratch/users/$USER/RecoDataDev/$jobname
    ln -nfs /pnfs/e1039/scratch/users/$USER/RecoDataDev data
else
    echo "Local mode."
    work=$dir_macros/scratch/$jobname
fi

mkdir -p $work
chmod -R 01755 $work

cd $dir_macros/..
tar -czvf $work/input.tar.gz  support RecoE1039Data.C
cd -

for run_num in $(cat $run_list) ; do  
  #location of the data in your area
  data_dir="/pnfs/e1039/scratch/users/$USER/DstRun"
  declare -a data_path_list=()
  if [ $dst_mode = 'single' ] ; then
      data_path_list=( $data_dir/$(printf 'run_%06d_spin.root' $run_num) )
  else # 'split'
      data_path_list=( $(find $data_dir -name $(printf 'run_%06d_spill_*_spin.root' $run_num) ) )
  fi

  for data_path in ${data_path_list[*]} ; do
    data_file=$(basename $data_path)
    job_name=${data_file%'.root'}
    mkdir -p $work/$job_name/log
    mkdir -p $work/$job_name/out
    chmod -R 01755 $work/$job_name
    
    rsync -av $dir_macros/gridrun_data.sh $work/$job_name/gridrun_data.sh

    if [ $do_sub == 1 ]; then
      CMD="/e906/app/software/script/jobsub_submit_spinquest.sh"
      CMD+=" --expected-lifetime='medium'" # medium=8h, short=3h, long=23h
      CMD+=" -L $work/$job_name/log/log.txt"
      CMD+=" -f $work/input.tar.gz"
      CMD+=" -d OUTPUT $work/$job_name/out"
      CMD+=" -f $data_path"
      CMD+=" file://$work/$job_name/gridrun_data.sh $nevents $run_num $data_file"
      echo "$CMD"
      unbuffer $CMD |& tee $work/$id/log_jobsub_submit.txt
      RET_SUB=${PIPESTATUS[0]}
      test $RET_SUB -ne 0 && exit $RET_SUB
    else
      mkdir -p $work/$job_name/input
      rsync -av $work/input.tar.gz $data_path  $work/$job_name/input
      cd $work/$job_name/
      $work/$job_name/gridrun_data.sh $nevents $run_num $data_file | tee $work/$job_name/log/log.txt
      cd -
    fi
  done
done

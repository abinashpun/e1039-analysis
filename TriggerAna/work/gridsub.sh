#!/bin/bash

dir_macros=$(dirname $(readlink -f $BASH_SOURCE))

jobname=$1
do_sub=$2
njobs=$3
nfile=$4
list=$5

echo "njobs=$njobs"
echo "nfile=$nfile"
echo "list=$list"

if [ $do_sub == 1 ]; then
    echo "Grid mode."
    dir_data=/pnfs/e1039/scratch/users/$USER/TriggerAna
    work=$dir_data/$jobname
    ln -nfs $dir_data data # for convenience
else
    echo "Local mode."
    work=$dir_macros/scratch/$jobname
fi

rm -rf         $work
mkdir -p       $work
chmod -R 01755 $work

cd $dir_macros
tar -czvf $work/input.tar.gz *.C *.cfg *.opts

for (( id=1; id<=$njobs; id++ ))
do  
  mkdir -p $work/$id/out
  chmod -R 01755 $work/$id
  cp -p $dir_macros/gridrun.sh $work/$id

  if [ $do_sub == 1 ]; then
    CMD="/e906/app/software/script/jobsub_submit_spinquest.sh"
    CMD+=" --expected-lifetime='medium'" # medium=8h, short=3h, long=23h
    CMD+=" -L $work/$id/log_gridrun.txt"
    CMD+=" -f $work/input.tar.gz"
    CMD+=" -d OUTPUT $work/$id/out"
    CMD+=" file://$work/$id/gridrun.sh $nfile $list $id"
    echo "$CMD"
    unbuffer $CMD |& tee $work/$id/log_jobsub_submit.txt
    RET_SUB=${PIPESTATUS[0]}
    test $RET_SUB -ne 0 && exit $RET_SUB
  else
    mkdir -p $work/$id/input
    cp -p $work/input.tar.gz $work/$id/input
    cd $work/$id
    $work/$id/gridrun.sh $nfile $list $id |& tee $work/$id/log_gridrun.txt
  fi
done

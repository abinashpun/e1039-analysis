#!/bin/bash

jobname=$1
do_sub=$2

njobs=$3
nevents=$4

gap=$5
target_pos=-300

if [ $do_sub == 1 ]; then
  echo "grid!"
else
  echo "local!"
fi

echo "njobs=$njobs"
echo "nevents=$nevents"

macros=/exp/seaquest/app/users/$USER/seaquest-analysis/HodoAccGap

sed "s/nevents=NAN/nevents=$nevents/"             $macros/gridrun.sh > $macros/gridrun_new.sh 
sed -i "s/gap=NAN/gap=$gap/"                      $macros/gridrun_new.sh
sed -i "s/target_pos=NAN/target_pos=$target_pos/" $macros/gridrun_new.sh
chmod +x $macros/gridrun_new.sh

if [ $do_sub == 1 ]; then
work=/pnfs/e1039/scratch/users/$USER/HodoAccGap/$jobname
else
work=$macros/scratch/$jobname
fi

mkdir -p $work
chmod -R 01755 $work

cd $macros
tar -czvf $work/input.tar.gz *.C *.cfg *.opts
cd -

for (( id=1; id<=$njobs; id++ ))
do  
  mkdir -p $work/$id/log
  mkdir -p $work/$id/out
  chmod -R 01755 $work/$id

  rsync -av $macros/gridrun_new.sh $work/$id/gridrun_new.sh

  CMD="/exp/seaquest/app/software/script/jobsub_submit_spinquest.sh"
  CMD+=" --expected-lifetime='medium'" # medium=8h, short=3h, long=23h
  CMD+=" -L $work/$id/log/log.txt"
  CMD+=" -f $work/input.tar.gz"
  CMD+=" -d OUTPUT $work/$id/out"
  CMD+=" file://$work/$id/gridrun_new.sh"

  if [ $do_sub == 1 ]; then
    echo $CMD
    $CMD
  else
    mkdir -p $work/$id/input
    rsync -av $work/input.tar.gz $work/$id/input
    cd $work/$id/
    $work/$id/gridrun_new.sh | tee $work/$id/log/log.txt
    cd -
  fi
done

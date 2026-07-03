#!/bin/bash

nevents=$1
run_num=$2

# simulate grid if not
is_grid=1
if [ -z ${CONDOR_DIR_INPUT+x} ];
  then
    CONDOR_DIR_INPUT=./input;
    echo "CONDOR_DIR_INPUT is initiallized as $CONDOR_DIR_INPUT"
    is_grid=0
  else
    echo "CONDOR_DIR_INPUT is set to '$CONDOR_DIR_INPUT'";
fi

if [ -z ${CONDOR_DIR_OUTPUT+x} ];
  then
    CONDOR_DIR_OUTPUT=./out;
    mkdir -p $CONDOR_DIR_OUTPUT
    echo "CONDOR_DIR_OUTPUT is initiallized as $CONDOR_DIR_OUTPUT"
  else
    echo "CONDOR_DIR_OUTPUT is set to '$CONDOR_DIR_OUTPUT'";
fi

# init grid log
echo "Hello, grid." | tee out.txt $CONDOR_DIR_OUTPUT/out.txt
pwd | tee -a out.txt $CONDOR_DIR_OUTPUT/out.txt

# setup input: macros; param folder; CODA data
tar -xzf $CONDOR_DIR_INPUT/input.tar.gz
tar -xzf $CONDOR_DIR_INPUT/para.tar.gz
ln -sf $CONDOR_DIR_INPUT/*.dat .
ls -lh | tee -a out.txt $CONDOR_DIR_OUTPUT/out.txt

# setup enviroment
#is_grid=1
if [ $is_grid == 1 ]; then
  echo "source /cvmfs/seaquest.opensciencegrid.org/seaquest/software/e1039/this-e1039.sh"
  source /cvmfs/seaquest.opensciencegrid.org/seaquest/software/e1039/this-e1039.sh
else
  echo "source /exp/seaquest/app/software/osg/software/e1039/this-e1039.sh"
  source /exp/seaquest/app/software/osg/software/e1039/this-e1039.sh
fi

# test enviroment
echo `which root`
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH"
echo "ldd $E1039_CORE/lib/libktracker.so"
ldd $E1039_CORE/lib/libktracker.so
echo "ldd $E1039_CORE/lib/libg4detectors.so"
ldd $E1039_CORE/lib/libg4detectors.so
echo "ldd $E1039_CORE/lib/libdecoder_maindaq.so"
ldd $E1039_CORE/lib/libdecoder_maindaq.so

# run
time root -b -q Fun4CODA.C\($nevents,$run_num\)

# mv output ROOT files to $CONDOR_DIR_OUTPUT
mv *.root $CONDOR_DIR_OUTPUT/

# signal success
echo "gridrun.sh finished!"

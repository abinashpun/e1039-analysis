DIR_NMR_DATA=/data2/e1039_data/target_data/NMR

alias list-nmr-event-files='ls $DIR_NMR_DATA/202*s.csv'

if [ -e /data2/e1039/this-e1039.sh ] ; then
    source /data2/e1039/this-e1039.sh
elif [ -e /exp/seaquest/app/software/osg/software/e1039/this-e1039.sh ] ; then
    source /exp/seaquest/app/software/osg/software/e1039/this-e1039.sh
else
    echo "Cannot find any this-e1039.sh.  You need set up ROOT by yourself."
fi

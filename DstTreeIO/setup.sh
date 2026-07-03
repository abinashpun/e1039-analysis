
if [ ${HOSTNAME:0:5} = 'e1039' ] ; then
    source /data2/e1039/this-e1039.sh 
elif [ ${HOSTNAME:0:13} = 'spinquestgpvm' ] ; then
    E1039_CORE_VERSION=pr.133
    source /exp/seaquest/app/software/osg/software/e1039/this-e1039.sh 
elif [ ${HOSTNAME:0:6} = 'udc-ba' ] ; then # UVA
    E1039_CORE_VERSION=pr.133
    source /project/ptgroup/spinquest/this-e1039.sh
fi

if [ -z $E1039_ROOT ] ; then
    echo "Cannot source 'this-e1039.sh' properly.  You could set up the E1039 environment by yourself."
fi

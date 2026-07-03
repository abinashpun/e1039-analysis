DIR_TOP=$(dirname $(readlink -f $BASH_SOURCE))
DIR_BUILD=$DIR_TOP/build
DIR_INST=$DIR_TOP/inst

## Choose one of e1039-core versions (or add your own one).
#E1039_CORE_DIR=/project/ptgroup/Jay/new-core/core-inst # Jay's modified ver.
E1039_CORE_VERSION=pr.133 # Latest system-wide ver.

source /project/ptgroup/spinquest/this-e1039.sh
#source /exp/seaquest/app/software/osg/software/e1039/this-e1039.sh
export LD_LIBRARY_PATH=$DIR_INST/lib:$LD_LIBRARY_PATH

function cmake-this {
    if [ -e $DIR_BUILD ] ; then
	echo "Clean up the build directory..."
	\rm -rf $DIR_BUILD
    fi
    if [ -e $DIR_INST ] ; then
	echo "Clean up the install directory..."
	\rm -rf $DIR_INST
    fi
    mkdir -p $DIR_BUILD
    echo "Run cmake..."
    ( cd $DIR_BUILD && cmake -DCMAKE_INSTALL_PREFIX=$DIR_INST $DIR_TOP/src )
    local RET=$?
    if [ $RET -eq 0 ] ; then
	echo "OK.  Run 'make-this' at any time/directory."
    else
	echo "NG.  Need a fix."
    fi
    return $RET
}

function make-this {
    if [ ! -e $DIR_BUILD/Makefile ] ; then
	echo "The build directory is not ready.  Complete 'cmake-this' first."
	return 1
    fi
    ( cd $DIR_BUILD && make install )
    return $?
}

##
## Main
##
echo "Two commands for source build:"
echo "  cmake-this"
echo "  make-this"
echo
echo "These commands can be executed in any directory."
echo "All built files are always created under $DIR_BUILD."
echo
echo "You need execute 'make-this'  when you modify source files."
echo "You need execute 'cmake-this' when you add/delete source files."

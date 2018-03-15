#!/bin/sh
cd ${0%/*} || exit 1                        # Run from this directory
. $WM_PROJECT_DIR/bin/tools/RunFunctions    # Tutorial run functions

# 1) Run meshing
# 2) Reconstruct
# 3) Test input zones and movement

#
# linkParallelCase srcDir dstDir
#
linkParallelCase()
{
    local src=$1
    local dst=$2
    shift 2

    if [ -e "$dst" ]
    then
        echo "Case already linked: remove case directory $dst prior to linking"
        return 1
    elif [ ! -d "$src" ]
    then
        echo "Error: no directory to link: $src"
        return 1
    fi

    echo "Linking $dst parallel case from $src"
    mkdir $dst

    # Copy system - may wish to change things
    for i in system 0
    do
        echo "    copy $i/"
        ( cd $dst && cp -r ../$src/$i . )
    done

    echo "    link constant/"
    ( cd $dst && ln -sf ../$src/constant . )

    echo "    link processor*/ with $# times: $@"
    for proc in $(cd $src && \ls -d processor*)
    do
        ( cd $dst && ln -sf ../$src/$proc . )
    done

    return 0
}


# Do steady-state case
(cd steady && ./Allrun.pre)

if notTest $@
then
    # Copy/link the steady-state case to movement
    linkParallelCase steady movement

    # Test movement
    \cp files/Allrun.movement movement/Allrun
    (cd movement && foamRunTutorials)
fi

#------------------------------------------------------------------------------

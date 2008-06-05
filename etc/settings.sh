#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software; you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation; either version 2 of the License, or (at your
#     option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM; if not, write to the Free Software Foundation,
#     Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
#
# Script
#     settings.sh
#
# Description
#     Startup file for OpenFOAM
#     Sourced from OpenFOAM-??/etc/bashrc
#
#------------------------------------------------------------------------------

if [ "$PS1" -a "$foamDotFile" ]; then
    if [ "$FOAM_VERBOSE" ]; then
        echo "Executing: $foamDotFile"
    fi
fi

AddPath()
{
    while [ $# -ge 1 ]
    do
        [ -d $1 ] || mkdir -p $1
        export PATH=$1:$PATH
        shift
    done
}

AddLib()
{
    while [ $# -ge 1 ]
    do
        [ -d $1 ] || mkdir -p $1
        export LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
        shift
    done
}


#- Add the system-specifc executables path to the path
export PATH=$WM_PROJECT_DIR/bin:$FOAM_INST_DIR/$WM_ARCH/bin:$PATH

#- Location of the jobControl directory
export FOAM_JOB_DIR=$FOAM_INST_DIR/jobControl

export WM_DIR=$WM_PROJECT_DIR/wmake
export WM_LINK_LANGUAGE=c++
export WM_OPTIONS=$WM_ARCH${WM_COMPILER}$WM_PRECISION_OPTION$WM_COMPILE_OPTION
export PATH=$WM_DIR:$PATH

export FOAM_SRC=$WM_PROJECT_DIR/src
export FOAM_LIB=$WM_PROJECT_DIR/lib
export FOAM_LIBBIN=$FOAM_LIB/$WM_OPTIONS
AddLib $FOAM_LIBBIN

export FOAM_APP=$WM_PROJECT_DIR/applications
export FOAM_APPBIN=$WM_PROJECT_DIR/applications/bin/$WM_OPTIONS
AddPath $FOAM_APPBIN

export FOAM_TUTORIALS=$WM_PROJECT_DIR/tutorials
export FOAM_UTILITIES=$FOAM_APP/utilities
export FOAM_SOLVERS=$FOAM_APP/solvers

export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/lib/$WM_OPTIONS
AddLib $FOAM_USER_LIBBIN

export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/applications/bin/$WM_OPTIONS
AddPath $FOAM_USER_APPBIN

export FOAM_RUN=$WM_PROJECT_USER_DIR/run


# Location of third-party software
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
thirdParty=$WM_PROJECT_INST_DIR/ThirdParty

# Compiler settings
# ~~~~~~~~~~~~~~~~~
WM_COMPILER_BIN=
WM_COMPILER_LIB=

# Select compiler installation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WM_COMPILER_INST = OpenFOAM | System
WM_COMPILER_INST=OpenFOAM

case "$WM_COMPILER_INST" in
OpenFOAM)
    case "$WM_COMPILER" in
    Gcc43)
        export WM_COMPILER_DIR=$thirdParty/gcc-4.3.0/platforms/$WM_ARCH$WM_COMPILER_ARCH
        ;;
    Gcc)
        export WM_COMPILER_DIR=$thirdParty/gcc-4.2.2/platforms/$WM_ARCH$WM_COMPILER_ARCH
        ;;
    esac

    # Check that the compiler directory can be found
    if [ ! -d "$WM_COMPILER_DIR" ]
    then
        echo
        echo "Warning in $WM_PROJECT_DIR/etc/settings.sh:"
        echo "    Cannot find $WM_COMPILER_DIR installation."
        echo "    Please install this compiler version or if you wish to use the system compiler,"
        echo "    change the WM_COMPILER_INST setting to 'System' in this file"
        echo
    fi

    WM_COMPILER_BIN=$WM_COMPILER_DIR/bin
    WM_COMPILER_LIB=$WM_COMPILER_DIR/lib$WM_COMPILER_LIB_ARCH:$WM_COMPILER_DIR/lib
    ;;
esac

if [ "$WM_COMPILER_BIN" != "" ]; then
    export PATH=$WM_COMPILER_BIN:$PATH
    export LD_LIBRARY_PATH=$WM_COMPILER_LIB:$LD_LIBRARY_PATH
fi


# Communications library
# ~~~~~~~~~~~~~~~~~~~~~~

unset MPI_ARCH_PATH

case "$WM_MPLIB" in
OPENMPI)
    mpi_version=openmpi-1.2.6
    export MPI_ARCH_PATH=$thirdParty/$mpi_version/platforms/$WM_OPTIONS

    # Tell OpenMPI where to find its install directory
    export OPAL_PREFIX=$MPI_ARCH_PATH

    AddLib  $MPI_ARCH_PATH/lib
    AddPath $MPI_ARCH_PATH/bin

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/$mpi_version
    unset mpi_version
    ;;

LAM)
    mpi_version=lam-7.1.4
    export MPI_ARCH_PATH=$thirdParty/$mpi_version/platforms/$WM_OPTIONS
    export LAMHOME=$thirdParty/$mpi_version
    # note: LAMHOME is deprecated, should probably point to MPI_ARCH_PATH too

    AddLib  $MPI_ARCH_PATH/lib
    AddPath $MPI_ARCH_PATH/bin

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/$mpi_version
    unset mpi_version
    ;;

MPICH)
    mpi_version=mpich-1.2.4
    export MPI_ARCH_PATH=$thirdParty/$mpi_version/platforms/$WM_OPTIONS
    export MPICH_ROOT=$MPI_ARCH_PATH

    AddLib  $MPI_ARCH_PATH/lib
    AddPath $MPI_ARCH_PATH/bin

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/$mpi_version
    unset mpi_version
    ;;

MPICH-GM)
    export MPI_ARCH_PATH=/opt/mpi
    export MPICH_PATH=$MPI_ARCH_PATH
    export MPICH_ROOT=$MPI_ARCH_PATH
    export GM_LIB_PATH=/opt/gm/lib64

    AddLib  $MPI_ARCH_PATH/lib
    AddLib  $GM_LIB_PATH
    AddPath $MPI_ARCH_PATH/bin

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/mpich-gm
    ;;

GAMMA)
    export MPI_ARCH_PATH=/usr
    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/gamma
    ;;

MPI)
    export MPI_ARCH_PATH=/opt/mpi
    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/mpi
    ;;

*)
    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/dummy
    ;;
esac

AddLib $FOAM_MPI_LIBBIN


# Set the MPI buffer size (used by all platforms except SGI MPI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export MPI_BUFFER_SIZE=20000000


# CGAL library if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~
[ -d "$CGAL_LIB_DIR" ] && AddLib $CGAL_LIB_DIR


# Switch on the hoard memory allocator if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if [ -f $FOAM_LIBBIN/libhoard.so ]; then
#    export LD_PRELOAD=$FOAM_LIBBIN/libhoard.so:$LD_PRELOAD
#fi

# -----------------------------------------------------------------------------

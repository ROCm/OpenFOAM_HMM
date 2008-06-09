#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
#    \\/     M anipulation  |
#-------------------------------------------------------------------------------
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
#     settings.csh
#
# Description
#     Startup file for OpenFOAM
#     Sourced from OpenFOAM-??/etc/cshrc
#
#------------------------------------------------------------------------------
if ($?prompt && $?foamDotFile) then
    if ($?FOAM_VERBOSE) then
        echo "Executing: $foamDotFile"
    endif
endif

alias AddPath 'set path=(\!* $path) ; if ( ! -d \!* ) mkdir -p \!*'
alias AddLib 'setenv LD_LIBRARY_PATH \!*\:${LD_LIBRARY_PATH} ; if ( ! -d \!* ) mkdir -p \!*'


#- Add the system-specific executables path to the path
set path=($WM_PROJECT_DIR/bin $FOAM_INST_DIR/$WM_ARCH/bin $path)

#- Location of the jobControl directory
setenv FOAM_JOB_DIR $FOAM_INST_DIR/jobControl

setenv WM_DIR $WM_PROJECT_DIR/wmake
setenv WM_LINK_LANGUAGE c++
setenv WM_OPTIONS $WM_ARCH$WM_COMPILER$WM_PRECISION_OPTION$WM_COMPILE_OPTION

set path=($WM_DIR $path)

setenv FOAM_SRC $WM_PROJECT_DIR/src
setenv FOAM_LIB $WM_PROJECT_DIR/lib
setenv FOAM_LIBBIN $FOAM_LIB/$WM_OPTIONS
AddLib $FOAM_LIBBIN

setenv FOAM_APP $WM_PROJECT_DIR/applications
setenv FOAM_APPBIN $WM_PROJECT_DIR/applications/bin/$WM_OPTIONS
AddPath $FOAM_APPBIN

setenv FOAM_TUTORIALS $WM_PROJECT_DIR/tutorials
setenv FOAM_UTILITIES $FOAM_APP/utilities
setenv FOAM_SOLVERS $FOAM_APP/solvers

setenv FOAM_USER_LIBBIN $WM_PROJECT_USER_DIR/lib/$WM_OPTIONS
AddLib $FOAM_USER_LIBBIN

setenv FOAM_USER_APPBIN $WM_PROJECT_USER_DIR/applications/bin/$WM_OPTIONS
AddPath $FOAM_USER_APPBIN

setenv FOAM_RUN $WM_PROJECT_USER_DIR/run


# Compiler settings
# ~~~~~~~~~~~~~~~~~
set WM_COMPILER_BIN=
set WM_COMPILER_LIB=

# Select compiler installation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WM_COMPILER_INST = OpenFOAM | System
set WM_COMPILER_INST=OpenFOAM


switch ("$WM_COMPILER_INST")
case OpenFOAM:
    switch ("$WM_COMPILER")
    case Gcc43:
        setenv WM_COMPILER_DIR $WM_THIRD_PARTY/gcc-4.3.0/platforms/$WM_ARCH$WM_COMPILER_ARCH
    breaksw
    case Gcc:
        setenv WM_COMPILER_DIR $WM_THIRD_PARTY/gcc-4.2.2/platforms/$WM_ARCH$WM_COMPILER_ARCH
    breaksw
    endsw

    # Check that the compiler directory can be found
    if ( ! -d "$WM_COMPILER_DIR" ) then
        echo
        echo "Warning in $WM_PROJECT_DIR/etc/settings.csh:"
        echo "    Cannot find $WM_COMPILER_DIR installation."
        echo "    Please install this compiler version or if you wish to use the system compiler,"
        echo "    change the WM_COMPILER_INST setting to 'System' in this file"
        echo
    endif

    set WM_COMPILER_BIN="$WM_COMPILER_DIR/bin"
    set WM_COMPILER_LIB=$WM_COMPILER_DIR/lib${WM_COMPILER_LIB_ARCH}:$WM_COMPILER_DIR/lib
    breaksw
endsw

if ($?WM_COMPILER_BIN) then
    set path=($WM_COMPILER_BIN $path)

    if ($?LD_LIBRARY_PATH) then
        setenv LD_LIBRARY_PATH ${WM_COMPILER_LIB}:${LD_LIBRARY_PATH}
    else
        setenv LD_LIBRARY_PATH ${WM_COMPILER_LIB}
    endif
endif


# Communications library
# ~~~~~~~~~~~~~~~~~~~~~~

unset MPI_ARCH_PATH

switch ("$WM_MPLIB")
case OPENMPI:
    set mpi_version=openmpi-1.2.6
    setenv MPI_HOME $WM_THIRD_PARTY/$mpi_version
    setenv MPI_ARCH_PATH $MPI_HOME/platforms/$WM_OPTIONS

    # Tell OpenMPI where to find its install directory
    setenv OPAL_PREFIX $MPI_ARCH_PATH

    AddLib  $MPI_ARCH_PATH/lib
    AddPath $MPI_ARCH_PATH/bin

    setenv FOAM_MPI_LIBBIN $FOAM_LIBBIN/$mpi_version
    unset mpi_version
    breaksw

case LAM:
    set mpi_version=lam-7.1.4
    setenv MPI_HOME $WM_THIRD_PARTY/$mpi_version
    setenv MPI_ARCH_PATH $MPI_HOME/platforms/$WM_OPTIONS
    setenv LAMHOME $WM_THIRD_PARTY/$mpi_version
    # note: LAMHOME is deprecated, should probably point to MPI_ARCH_PATH too

    AddLib  $MPI_ARCH_PATH/lib
    AddPath $MPI_ARCH_PATH/bin

    setenv FOAM_MPI_LIBBIN $FOAM_LIBBIN/$mpi_version
    unset mpi_version
    breaksw

case MPICH:
    set mpi_version=mpich-1.2.4
    setenv MPI_HOME $WM_THIRD_PARTY/$mpi_version
    setenv MPI_ARCH_PATH $MPI_HOME/platforms/$WM_OPTIONS
    setenv MPICH_ROOT $MPI_ARCH_PATH

    AddLib  $MPI_ARCH_PATH/lib
    AddPath $MPI_ARCH_PATH/bin

    setenv FOAM_MPI_LIBBIN $FOAM_LIBBIN/$mpi_version
    unset mpi_version
    breaksw

case MPICH-GM:
    setenv MPI_ARCH_PATH /opt/mpi
    setenv MPICH_PATH $MPI_ARCH_PATH
    setenv MPICH_ROOT $MPI_ARCH_PATH
    setenv GM_LIB_PATH /opt/gm/lib64

    AddLib  $MPI_ARCH_PATH/lib
    AddLib  $GM_LIB_PATH
    AddPath $MPI_ARCH_PATH/bin

    setenv FOAM_MPI_LIBBIN $FOAM_LIBBIN/mpich-gm
    breaksw

case GAMMA:
    setenv MPI_ARCH_PATH /usr
    setenv FOAM_MPI_LIBBIN $FOAM_LIBBIN/gamma
    breaksw

case MPI:
    setenv MPI_ARCH_PATH /opt/mpi
    setenv FOAM_MPI_LIBBIN $FOAM_LIBBIN/mpi
    breaksw

default:
    setenv FOAM_MPI_LIBBIN $FOAM_LIBBIN/dummy
    breaksw
endsw

AddLib $FOAM_MPI_LIBBIN


# Set the MPI buffer size (used by all platforms except SGI MPI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setenv MPI_BUFFER_SIZE 20000000


# CGAL library if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~
if ( $?CGAL_LIB_DIR ) then
    AddLib $CGAL_LIB_DIR
endif


# Switch on the hoard memory allocator if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if ( -f $FOAM_LIBBIN/libhoard.so ) then
#    setenv LD_PRELOAD $FOAM_LIBBIN/libhoard.so:${LD_PRELOAD}
#endif

# -----------------------------------------------------------------------------

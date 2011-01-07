#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     etc/settings.csh
#
# Description
#     Startup file for OpenFOAM
#     Sourced from OpenFOAM-<VERSION>/etc/cshrc
#
#------------------------------------------------------------------------------

# prefix to PATH
alias _foamAddPath 'setenv PATH \!*\:${PATH}'
# prefix to LD_LIBRARY_PATH
alias _foamAddLib 'setenv LD_LIBRARY_PATH \!*\:${LD_LIBRARY_PATH}'
# prefix to MANPATH
alias _foamAddMan 'setenv MANPATH \!*\:${MANPATH}'

#------------------------------------------------------------------------------

# location of the jobControl directory
setenv FOAM_JOB_DIR $WM_PROJECT_INST_DIR/jobControl

# wmake configuration
setenv WM_DIR $WM_PROJECT_DIR/wmake
setenv WM_LINK_LANGUAGE c++
setenv WM_OPTIONS $WM_ARCH$WM_COMPILER$WM_PRECISION_OPTION$WM_COMPILE_OPTION

# base executables/libraries
setenv FOAM_APPBIN $WM_PROJECT_DIR/bin/$WM_OPTIONS
setenv FOAM_LIBBIN $WM_PROJECT_DIR/lib/$WM_OPTIONS

# external (ThirdParty) libraries
setenv FOAM_EXT_LIBBIN $WM_THIRD_PARTY_DIR/lib/$WM_OPTIONS

# shared site executables/libraries
# similar naming convention as ~OpenFOAM expansion
setenv FOAM_SITE_APPBIN $WM_PROJECT_INST_DIR/site/$WM_PROJECT_VERSION/bin/$WM_OPTIONS
setenv FOAM_SITE_LIBBIN $WM_PROJECT_INST_DIR/site/$WM_PROJECT_VERSION/lib/$WM_OPTIONS

# user executables/libraries
setenv FOAM_USER_APPBIN $WM_PROJECT_USER_DIR/bin/$WM_OPTIONS
setenv FOAM_USER_LIBBIN $WM_PROJECT_USER_DIR/lib/$WM_OPTIONS

# convenience
setenv FOAM_APP $WM_PROJECT_DIR/applications
setenv FOAM_LIB $WM_PROJECT_DIR/lib
setenv FOAM_SRC $WM_PROJECT_DIR/src
setenv FOAM_TUTORIALS $WM_PROJECT_DIR/tutorials
setenv FOAM_UTILITIES $FOAM_APP/utilities
setenv FOAM_SOLVERS $FOAM_APP/solvers
setenv FOAM_RUN $WM_PROJECT_USER_DIR/run

# add wmake to the path - not required for runtime only environment
if ( -d "${WM_DIR}" ) setenv PATH ${WM_DIR}:${PATH}
# add OpenFOAM scripts to the path
setenv PATH ${WM_PROJECT_DIR}/bin:${PATH}

_foamAddPath ${FOAM_USER_APPBIN}:${FOAM_SITE_APPBIN}:${FOAM_APPBIN}
# Make sure to pick up dummy versions of external libraries last
_foamAddLib  ${FOAM_USER_LIBBIN}:${FOAM_SITE_LIBBIN}:${FOAM_LIBBIN}:${FOAM_EXT_LIBBIN}:${FOAM_LIBBIN}/dummy

# Compiler settings
# ~~~~~~~~~~~~~~~~~
unset gcc_version gmp_version mpfr_version mpc_version
unsetenv MPFR_ARCH_PATH


# Location of compiler installation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( ! $?foamCompiler ) then
then
    foamCompiler=system
    echo "Warning in $WM_PROJECT_DIR/etc/settings.csh:"
    echo "    foamCompiler not set, using '$foamCompiler'"
endif

switch ("$foamCompiler")
case OpenFOAM:
case ThirdParty:
    switch ("$WM_COMPILER")
    case Gcc:
    case Gcc++0x:
        set gcc_version=gcc-4.4.3
        set gmp_version=gmp-5.0.1
        set mpfr_version=mpfr-2.4.2
        breaksw
    case Gcc45:
    case Gcc45++0x:
        set gcc_version=gcc-4.5.2
        set gmp_version=gmp-5.0.1
        set mpfr_version=mpfr-2.4.2
        set mpc_version=mpc-0.8.1
        breaksw
    case Gcc44:
    case Gcc44++0x:
        set gcc_version=gcc-4.4.3
        set gmp_version=gmp-5.0.1
        set mpfr_version=mpfr-2.4.2
        breaksw
    case Gcc43:
        set gcc_version=gcc-4.3.3
        set gmp_version=gmp-4.2.4
        set mpfr_version=mpfr-2.4.1
        breaksw
    case Clang:
        # using clang - not gcc
        setenv WM_CC 'clang'
        setenv WM_CXX 'clang++'
        #set clang_version=llvm-2.8
        set clang_version=llvm-svn
        breaksw
    default:
        echo
        echo "Warning in $WM_PROJECT_DIR/etc/settings.csh:"
        echo "    Unknown OpenFOAM compiler type '$WM_COMPILER'"
        echo "    Please check your settings"
        echo
        breaksw
    endsw

    if ( $?gcc_version ) then
        set gccDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$gcc_version
        set gmpDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$gmp_version
        set mpfrDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$mpfr_version
        if ( $?mpc_version ) then
            set mpcDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$mpc_version
        endif

        # Check that the compiler directory can be found
        if ( ! -d "$gccDir" ) then
            echo
            echo "Warning in $WM_PROJECT_DIR/etc/settings.csh:"
            echo "    Cannot find $gccDir installation."
            echo "    Please install this compiler version or if you wish to use the system compiler,"
            echo "    change the 'foamCompiler' setting to 'system'"
            echo
        endif

        _foamAddMan     $gccDir/man
        _foamAddPath    $gccDir/bin

        # 64-bit needs lib64, but 32-bit needs lib (not lib32)
        if ($WM_ARCH_OPTION == 64) then
            _foamAddLib     $gccDir/lib$WM_COMPILER_LIB_ARCH
        else
            _foamAddLib     $gccDir/lib
        endif

        # add in gmp/mpfr libraries
        _foamAddLib     $gmpDir/lib
        _foamAddLib     $mpfrDir/lib

        # add in mpc libraries (not need for older gcc)
        if ( $?mpc_version ) then
            _foamAddLib     $mpcDir/lib
        endif

        # used by boost/CGAL:
        setenv MPFR_ARCH_PATH $mpfrDir
    endif
    unset gcc_version gccDir
    unset gmp_version gmpDir  mpfr_version mpfrDir  mpc_version mpcDir

    if ( $?clang_version ) then
        set clangDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$clang_version

        # Check that the compiler directory can be found
        if ( ! -d "$clangDir" ) then
            echo
            echo "Warning in $WM_PROJECT_DIR/etc/settings.csh:"
            echo "    Cannot find $clangDir installation."
            echo "    Please install this compiler version or if you wish to use the system compiler,"
            echo "    change the 'foamCompiler' setting to 'system'"
            echo
        endif

        _foamAddMan     $clangDir/man
        _foamAddPath    $clangDir/bin
    endif
    unset clang_version clangDir

    breaksw

case system:
    # okay, use system compiler
    breaksw

default:
    echo "Warn: foamCompiler='$foamCompiler' is unsupported"
    echo "   treating as 'system' instead"
    breaksw
endsw


#
# add c++0x flags for external programs
#
if ( $?WM_CXXFLAGS ) then
    switch ("$WM_COMPILER")
    case Gcc*++0x:
        setenv WM_CXXFLAGS "$WM_CXXFLAGS -std=c++0x"
        breaksw
    endsw
endif


# boost and CGAL
# ~~~~~~~~~~~~~~

set boost_version=boost_1_42_0
set cgal_version=CGAL-3.5.1

setenv BOOST_ARCH_PATH $WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$boost_version
setenv CGAL_ARCH_PATH  $WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$cgal_version

# enabled if CGAL is available
if ( -d "$CGAL_ARCH_PATH" ) then
    if ( -d "$BOOST_ARCH_PATH" ) then
        _foamAddLib $BOOST_ARCH_PATH/lib
    else
        unsetenv BOOST_ARCH_PATH
    endif
    _foamAddLib $CGAL_ARCH_PATH/lib
else
    unsetenv BOOST_ARCH_PATH CGAL_ARCH_PATH MPFR_ARCH_PATH
endif

unset boost_version cgal_version


# Communications library
# ~~~~~~~~~~~~~~~~~~~~~~

unsetenv MPI_ARCH_PATH MPI_HOME

switch ("$WM_MPLIB")
case OPENMPI:
    #set mpi_version=openmpi-1.4.1
    set mpi_version=openmpi-1.5.1
    setenv MPI_ARCH_PATH $WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$mpi_version

    # Tell OpenMPI where to find its install directory
    setenv OPAL_PREFIX $MPI_ARCH_PATH

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib
    _foamAddMan     $MPI_ARCH_PATH/man

    setenv FOAM_MPI_LIBBIN $FOAM_EXT_LIBBIN/$mpi_version
    unset mpi_version
    breaksw

case SYSTEMOPENMPI:
    # Use the system installed openmpi, get library directory via mpicc
    set mpi_version=openmpi-system

    # Set compilation flags here instead of in wmake/rules/../mplibSYSTEMOPENMPI
    setenv PINC "`mpicc --showme:compile`"
    setenv PLIBS "`mpicc --showme:link`"
    set libDir=`echo "$PLIBS" | sed -e 's/.*-L\([^ ]*\).*/\1/'`

    if ($?FOAM_VERBOSE && $?prompt) then
        echo "Using system installed MPI:"
        echo "    compile flags : $PINC"
        echo "    link flags    : $PLIBS"
        echo "    libmpi dir    : $libDir"
    endif

    _foamAddLib     $libDir

    setenv FOAM_MPI_LIBBIN $FOAM_EXT_LIBBIN/$mpi_version
    unset mpi_version libDir
    breaksw

case MPICH:
    set mpi_version=mpich2-1.1.1p1
    setenv MPI_HOME $WM_THIRD_PARTY_DIR/$mpi_version
    setenv MPI_ARCH_PATH $WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$mpi_version

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib
    _foamAddMan     $MPI_ARCH_PATH/share/man

    setenv FOAM_MPI_LIBBIN $FOAM_EXT_LIBBIN/$mpi_version
    unset mpi_version
    breaksw

case MPICH-GM:
    setenv MPI_ARCH_PATH /opt/mpi
    setenv MPICH_PATH $MPI_ARCH_PATH
    setenv GM_LIB_PATH /opt/gm/lib64

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib
    _foamAddLib     $GM_LIB_PATH

    setenv FOAM_MPI_LIBBIN $FOAM_EXT_LIBBIN/mpich-gm
    breaksw

case HPMPI:
    setenv MPI_HOME /opt/hpmpi
    setenv MPI_ARCH_PATH $MPI_HOME

    _foamAddPath $MPI_ARCH_PATH/bin

    switch (`uname -m`)
    case i686:
        _foamAddLib $MPI_ARCH_PATH/lib/linux_ia32
        breaksw
    case x86_64:
        _foamAddLib $MPI_ARCH_PATH/lib/linux_amd64
        breaksw
    case ia64:
        _foamAddLib $MPI_ARCH_PATH/lib/linux_ia64
        breaksw
    default:
        echo Unknown processor type `uname -m` for Linux
        breaksw
    endsw

    setenv FOAM_MPI_LIBBIN $FOAM_EXT_LIBBIN/hpmpi
    breaksw

case GAMMA:
    setenv MPI_ARCH_PATH /usr
    setenv FOAM_MPI_LIBBIN $FOAM_EXT_LIBBIN/gamma
    breaksw

case MPI:
    setenv MPI_ARCH_PATH /opt/mpi
    setenv FOAM_MPI_LIBBIN $FOAM_EXT_LIBBIN/mpi
    breaksw

case FJMPI:
    setenv MPI_ARCH_PATH /opt/FJSVmpi2
    setenv FOAM_MPI_LIBBIN $FOAM_EXT_LIBBIN/mpi
    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib/sparcv9
    _foamAddLib     /opt/FSUNf90/lib/sparcv9
    _foamAddLib     /opt/FJSVpnidt/lib
    breaksw

case QSMPI:
    setenv MPI_ARCH_PATH /usr/lib/mpi
    setenv FOAM_MPI_LIBBIN FOAM_LIBBIN/qsmpi

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib

    breaksw

default:
    setenv FOAM_MPI_LIBBIN $FOAM_LIBBIN/dummy
    breaksw
endsw

_foamAddLib $FOAM_MPI_LIBBIN


# Set the minimum MPI buffer size (used by all platforms except SGI MPI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( ! $?minBufferSize ) set minBufferSize=20000000

if ( $?MPI_BUFFER_SIZE ) then
    if ( $MPI_BUFFER_SIZE < $minBufferSize ) then
        setenv MPI_BUFFER_SIZE $minBufferSize
    endif
else
    setenv MPI_BUFFER_SIZE $minBufferSize
endif


# Enable the hoard memory allocator if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if ( -f $FOAM_EXT_LIBBIN/libhoard.so ) then
#    setenv LD_PRELOAD $FOAM_EXT_LIBBIN/libhoard.so:$LD_PRELOAD
#endif


# cleanup environment:
# ~~~~~~~~~~~~~~~~~~~~
unalias _foamAddPath _foamAddLib _foamAddMan
unset foamCompiler minBufferSize

# ----------------------------------------------------------------- end-of-file

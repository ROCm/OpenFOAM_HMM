#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, licensed under GNU General Public License
#     <http://www.gnu.org/licenses/>.
#
# File
#     config.csh/example/prefs.csh
#
# Description
#     Preset variables for the OpenFOAM configuration - C-Shell shell syntax.
#
#     The prefs.csh file will be sourced by the OpenFOAM etc/cshrc when it is
#     found by foamEtcFile.
#
# See also
#     'foamEtcFile -help' or 'foamEtcFile -list' for information about the
#     paths searched
#
#------------------------------------------------------------------------------

#- Compiler location:
setenv WM_COMPILER_TYPE ThirdParty

#- Compiler:
setenv WM_COMPILER Clang

#- MPI implementation:
setenv WM_MPLIB SYSTEMOPENMPI

#------------------------------------------------------------------------------

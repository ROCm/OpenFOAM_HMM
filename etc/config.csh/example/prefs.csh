#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | www.openfoam.com
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2011-2016 OpenFOAM Foundation
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, licensed under GNU General Public License
#     <http://www.gnu.org/licenses/>.
#
# File
#     config.csh/example/prefs.csh
#     - sourced by OpenFOAM-*/etc/cshrc
#
# Description
#     Example of preset variables for the OpenFOAM configuration (C-Shell shell)
#
# See also
#     'foamEtcFile -help' or 'foamEtcFile -list' for information about the
#     paths searched
#
#------------------------------------------------------------------------------

setenv WM_COMPILER_TYPE ThirdParty
setenv WM_COMPILER Clang
setenv WM_MPLIB SYSTEMOPENMPI

#------------------------------------------------------------------------------

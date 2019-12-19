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
#     config.sh/example/prefs.sh
#     - sourced by OpenFOAM-*/etc/bashrc
#
# Description
#     Example of preset variables for the OpenFOAM configuration (POSIX shell)
#
# See also
#     'foamEtcFile -help' or 'foamEtcFile -list' for information about the
#     paths searched
#
#------------------------------------------------------------------------------

export WM_COMPILER_TYPE=ThirdParty
export WM_COMPILER=Clang
export WM_MPLIB=SYSTEMOPENMPI

#------------------------------------------------------------------------------

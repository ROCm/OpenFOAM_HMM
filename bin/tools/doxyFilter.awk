#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
#    \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
# Script
#     doxyFilter.awk
#
# Description
#     Converts cocoon style sentinel strings into doxygen style strings
#
#     Assumes comment strings are formatted as follows
#         //- General description
#         //  Detailed information
#         //  and even more information
#     or
#         //- General description
#         //- that spans several
#         //- lines
#         //  Detailed information
#         //  and even more information
#
#     This should be re-formatted as the following
#         /*! \brief General description
#          that spans several
#          lines
#         */
#         /*! Detailed information
#         and even more information
#         */
#
#     The intermediate "/*! ... */" block is left-justified to handle
#     possible verbatim text
#
#------------------------------------------------------------------------------

# States: 0=normal, 1=brief, 2=details
BEGIN {
    state = 0
}

/^ *\/\/-/ {
    if (state == 0)
    {
        # Changed from normal to brief (start of comment block)
        printf "/*! \\brief"
        state = 1
    }

    if (state == 1)
    {
        # Within brief: strip leading
        if (!sub(/^ *\/\/-  /, ""))
        {
            sub(/^ *\/\/-/, "")
        }
    }

    print
    next
}

/^ *\/\// {
    if (state == 1)
    {
        # Change from brief to details
        printf "*/\n"
        printf "/*! "
        state = 2
    }

    if (state == 2)
    {
        # Within details: strip leading
        if (!sub(/^ *\/\/  /, ""))
        {
            sub(/^ *\/\//, "")
        }
    }

    print
    next
}

{
    # End comment filtering
    if (state)
    {
        printf "*/\n"
    }
    state = 0
    print
    next
}

#------------------------------------------------------------------------------

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Namespace
    Foam::stringOps

Description
    Collection of static functions to do various simple string-related
    operations

SourceFiles
    stringOps.C

\*---------------------------------------------------------------------------*/
#ifndef stringOps_H
#define stringOps_H

#include "string.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Namespace stringOps Declaration
\*---------------------------------------------------------------------------*/

namespace stringOps
{
    //- Expand occurences of variables according to the mapping
    //  Expansion includes:
    //  -# variables
    //    - "$VAR", "${VAR}"
    //
    //  \note the leading sigil can be changed to avoid conflicts with other
    //  string expansions
    string expand
    (
        const string&,
        const HashTable<string, word, string::hash>& mapping,
        const char sigil = '$'
    );


    //- Inplace expand occurences of variables according to the mapping
    //  Expansion includes:
    //  -# variables
    //    - "$VAR", "${VAR}"
    //
    //  \note the leading sigil can be changed to avoid conflicts with other
    //  string expansions
    string& inplaceExpand
    (
        string&,
        const HashTable<string, word, string::hash>& mapping,
        const char sigil = '$'
    );


    //- Expand initial tildes and all occurences of environment variables
    //  Expansion includes:
    //  -# environment variables
    //    - "$VAR", "${VAR}"
    //  -# current directory
    //    - leading "./" : the current directory
    //  -# tilde expansion
    //    - leading "~/" : home directory
    //    - leading "~user" : home directory for specified user
    //    - leading "~OpenFOAM" : site/user OpenFOAM configuration directory
    //
    //  \sa
    //  Foam::findEtcFile
    string expandEnv
    (
        const string&,
        const bool recurse=false,
        const bool allowEmptyVar = false
    );


    //- Expand initial tildes and all occurences of environment variables
    //  Expansion includes:
    //  -# environment variables
    //    - "$VAR", "${VAR}"
    //  -# current directory
    //    - leading "./" : the current directory
    //  -# tilde expansion
    //    - leading "~/" : home directory
    //    - leading "~user" : home directory for specified user
    //    - leading "~OpenFOAM" : site/user OpenFOAM configuration directory
    //
    //  \sa
    //  Foam::findEtcFile
    string& inplaceExpandEnv
    (
        string&,
        const bool recurse=false,
        const bool allowEmptyVar = false
    );


    //- Return string trimmed of leading whitespace
    string trimLeft(const string&);

    //- Trim leading whitespace inplace
    string& inplaceTrimLeft(string&);

    //- Return string trimmed of trailing whitespace
    //  NOT IMPLEMENTED
    string trimRight(const string&);

    //- Trim trailing whitespace inplace
    //  NOT IMPLEMENTED
    string& inplaceTrimRight(string&);

    //- Return string trimmed of leading and trailing whitespace
    //  NOT IMPLEMENTED
    string trim(const string&);

    //- Trim leading and trailing whitespace inplace
    //  NOT IMPLEMENTED
    string& inplaceTrim(string&);



} // End namespace stringOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

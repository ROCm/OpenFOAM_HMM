/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include <cctype>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::ensight::VarName::VarName(const char* s)
:
    word(s, false)
{
    stripInvalid();
}


inline Foam::ensight::VarName::VarName(const std::string& s)
:
    word(s, false)
{
    stripInvalid();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::ensight::VarName::valid(char c)
{
    return
    (
        word::valid(c) // includes space, quotes, /\;{}
     && c != '!'
     && c != '#'
     && c != '%'
     && c != '('
     && c != ')'
     && c != '*'
     && c != '+'
     && c != ','
     && c != '-'
     && c != '.'
     && c != '@'
     && c != '['
     && c != ']'
     && c != '^'
    );
}


inline void Foam::ensight::VarName::stripInvalid()
{
    string::stripInvalid<VarName>(*this);

    if (empty())
    {
        FatalErrorInFunction
            << "ensight::VarName empty after stripping" << nl
            << exit(FatalError);
    }

    // Prefix with '_' to avoid starting with leading digits
    std::string::iterator iter = begin();
    if (isdigit(*iter))
    {
        insert(iter, '_');
    }
}


// ************************************************************************* //

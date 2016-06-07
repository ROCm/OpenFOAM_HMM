/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include <cctype>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::ensight::FileName::FileName(const FileName& fn)
:
    fileName(fn, false)
{}


inline Foam::ensight::FileName::FileName(const char* s)
:
    fileName(s, false)
{
    stripInvalid();
}


inline Foam::ensight::FileName::FileName(const std::string& s)
:
    fileName(s, false)
{
    stripInvalid();
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::ensight::FileName::stripInvalid()
{
    string::stripInvalid<FileName>(*this);

    removeRepeated('/');
    removeTrailing('/');

    if (empty())
    {
        FatalErrorInFunction
            << "ensight::FileName empty after stripping" << nl
            << exit(FatalError);
    }
}


inline bool Foam::ensight::FileName::valid(char c)
{
    return
    (
        fileName::valid(c)  // includes space, quotes
     && c != '*'   // wild-card
     && c != '%'   // structured block continuation
    );
}


// ************************************************************************* //

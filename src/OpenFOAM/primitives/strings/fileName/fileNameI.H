/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include <iostream>  // for std::cerr

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::fileName::stripInvalid()
{
    // Skip stripping unless debug is active (to avoid costly operations)
    if (debug && string::stripInvalid<fileName>(*this))
    {
        std::cerr
            << "fileName::stripInvalid() called for invalid fileName "
            << this->c_str() << std::endl;

        if (debug > 1)
        {
            std::cerr
                << "    For debug level (= " << debug
                << ") > 1 this is considered fatal" << std::endl;
            std::abort();
        }

        removeRepeated('/');
        removeTrailing('/');
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::fileName::fileName()
:
    string()
{}


inline Foam::fileName::fileName(const fileName& fn)
:
    string(fn)
{}


inline Foam::fileName::fileName(const word& w)
:
    string(w)
{}


inline Foam::fileName::fileName(const string& s, const bool doStripInvalid)
:
    string(s)
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


inline Foam::fileName::fileName(const std::string& s, const bool doStripInvalid)
:
    string(s)
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


inline Foam::fileName::fileName(const char* s, const bool doStripInvalid)
:
    string(s)
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::fileName::valid(char c)
{
    return
    (
        !isspace(c)
     && c != '"'   // string quote
     && c != '\''  // string quote
    );
}


inline bool Foam::fileName::isAbsolute(const std::string& str)
{
    return !str.empty() && str[0] == '/';
}


inline bool Foam::fileName::isAbsolute() const
{
    return !empty() && operator[](0) == '/';
}


inline bool Foam::fileName::hasExt() const
{
    return string::hasExt();
}


inline bool Foam::fileName::removeExt()
{
    return string::removeExt();
}


// ************************************************************************* //

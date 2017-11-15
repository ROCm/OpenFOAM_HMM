/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include <cctype>
#include <iostream>  // for std::cerr

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::word::stripInvalid()
{
    // Skip stripping unless debug is active (to avoid costly operations)
    if (debug && string::stripInvalid<word>(*this))
    {
        std::cerr
            << "word::stripInvalid() called for word "
            << this->c_str() << std::endl;

        if (debug > 1)
        {
            std::cerr
                << "    For debug level (= " << debug
                << ") > 1 this is considered fatal" << std::endl;
            std::abort();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::word::word()
:
    string()
{}


inline Foam::word::word(const word& w)
:
    string(w)
{}


inline Foam::word::word(const string& s, const bool doStripInvalid)
:
    string(s)
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


inline Foam::word::word(const std::string& s, const bool doStripInvalid)
:
    string(s)
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


inline Foam::word::word(const char* s, const bool doStripInvalid)
:
    string(s)
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


inline Foam::word::word
(
    const char* s,
    const size_type len,
    const bool doStripInvalid
)
:
    string(s, len)
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


inline Foam::word::word(word&& w)
:
    string(std::move(w))
{}


inline Foam::word::word(string&& s, const bool doStripInvalid)
:
    string(std::move(s))
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


inline Foam::word::word(std::string&& s, const bool doStripInvalid)
:
    string(std::move(s))
{
    if (doStripInvalid)
    {
        stripInvalid();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::word::valid(char c)
{
    return
    (
        !isspace(c)
     && c != '"'   // string quote
     && c != '\''  // string quote
     && c != '/'   // path separator
     && c != ';'   // end statement
     && c != '{'   // beg subdict
     && c != '}'   // end subdict
    );
}


inline bool Foam::word::hasExt() const
{
    return string::hasExt();
}


inline bool Foam::word::removeExt()
{
    return string::removeExt();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::word::operator=(const word& w)
{
    string::operator=(w);
}


inline void Foam::word::operator=(word&& w)
{
    string::operator=(std::move(w));
}


inline void Foam::word::operator=(const string& s)
{
    string::operator=(s);
    stripInvalid();
}


inline void Foam::word::operator=(string&& s)
{
    string::operator=(std::move(s));
    stripInvalid();
}


inline void Foam::word::operator=(const std::string& s)
{
    string::operator=(s);
    stripInvalid();
}


inline void Foam::word::operator=(std::string&& s)
{
    string::operator=(std::move(s));
    stripInvalid();
}


inline void Foam::word::operator=(const char* s)
{
    string::operator=(s);
    stripInvalid();
}


// ************************************************************************* //

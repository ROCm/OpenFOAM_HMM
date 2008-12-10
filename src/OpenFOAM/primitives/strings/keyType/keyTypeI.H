/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::keyType::keyType()
:
    word(),
    isWildCard_(false)
{}


inline Foam::keyType::keyType(const keyType& s)
:
    word(s, false),
    isWildCard_(s.isPattern())
{}


inline Foam::keyType::keyType(const word& s)
:
    word(s, false),
    isWildCard_(false)
{}


// Construct as copy of string. Expect it to be regular expression
inline Foam::keyType::keyType(const string& s)
:
    word(s, false),
    isWildCard_(true)
{}


// Construct as copy of character array
inline Foam::keyType::keyType(const char* s)
:
    word(s, false),
    isWildCard_(false)
{}


//- Construct as copy of std::string
inline Foam::keyType::keyType
(
    const std::string& s,
    const bool isPattern
)
:
    word(s, false),
    isWildCard_(isPattern)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::keyType::valid(char c)
{
    return c != '"';
}


bool Foam::keyType::isPattern() const
{
    return isWildCard_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::keyType::operator=(const keyType& s)
{
    // Bypass checking
    string::operator=(s);
    isWildCard_ = s.isPattern();
}


inline void Foam::keyType::operator=(const word& s)
{
    word::operator=(s);
    isWildCard_ = false;
}


inline void Foam::keyType::operator=(const string& s)
{
    // Bypass checking
    string::operator=(s);
    isWildCard_ = true;
}


inline void Foam::keyType::operator=(const char* s)
{
    // Bypass checking
    string::operator=(s);
    isWildCard_ = false;
}


// ************************************************************************* //

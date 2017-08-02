/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline bool Foam::keyType::valid(char c)
{
    return
    (
        !isspace(c)
     && c != '"'   // string quote
     && c != '\''  // string quote
     && c != '/'   // path separator
     && c != ';'   // end statement
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::keyType::keyType()
:
    word(),
    isPattern_(false)
{}


inline Foam::keyType::keyType(const keyType& s)
:
    word(s, false),
    isPattern_(s.isPattern())
{}


inline Foam::keyType::keyType(const word& s)
:
    word(s, false),
    isPattern_(false)
{}


inline Foam::keyType::keyType(const string& s)
:
    word(s, false),
    isPattern_(true)
{}


inline Foam::keyType::keyType(const char* s)
:
    word(s, false),
    isPattern_(false)
{}


inline Foam::keyType::keyType
(
    const std::string& s,
    const bool isPattern
)
:
    word(s, false),
    isPattern_(isPattern)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::keyType::isPattern() const
{
    return isPattern_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::keyType::operator()(const std::string& text) const
{
    return match(text);  // Use smart match
}


inline void Foam::keyType::operator=(const keyType& s)
{
    string::operator=(s); // Bypass checking
    isPattern_ = s.isPattern_;
}


inline void Foam::keyType::operator=(const word& s)
{
    string::operator=(s); // Bypass checking
    isPattern_ = false;
}


inline void Foam::keyType::operator=(const string& s)
{
    string::operator=(s); // Bypass checking
    isPattern_ = true;
}


inline void Foam::keyType::operator=(const char* s)
{
    string::operator=(s); // Bypass checking
    isPattern_ = false;
}


// ************************************************************************* //

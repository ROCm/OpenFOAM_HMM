/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline bool Foam::regExp::meta(const char c)
{
    return
    (
        (c == '.')                           // any character
     || (c == '*' || c == '+' || c == '?')   // quantifiers
     || (c == '(' || c == ')' || c == '|')   // grouping/branching
     || (c == '[' || c == ']')               // range
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline bool Foam::regExp::empty() const
{
    return !preg_;
}


inline bool Foam::regExp::exists() const
{
    return preg_ ? true : false;
}


inline unsigned Foam::regExp::ngroups() const
{
    return preg_ ? preg_->re_nsub : 0;
}


inline bool Foam::regExp::search(const std::string& text) const
{
    return std::string::npos != find(text);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

inline void Foam::regExp::operator=(const char* pattern)
{
    set(pattern);
}


inline void Foam::regExp::operator=(const std::string& pattern)
{
    set(pattern);
}


// ************************************************************************* //

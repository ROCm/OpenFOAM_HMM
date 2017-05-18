/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include <iostream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::string::string()
{}


inline Foam::string::string(const std::string& str)
:
    std::string(str)
{}


// Copy character array
inline Foam::string::string(const char* str)
:
    std::string(str)
{}


// Construct from a given number of characters in a character array
inline Foam::string::string(const char* str, const size_type len)
:
    std::string(str, len)
{}


// Construct from a single character
inline Foam::string::string(const char c)
:
    std::string(1, c)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class String>
inline bool Foam::string::valid(const std::string& str)
{
    for (const_iterator iter = str.begin(); iter != str.end(); ++iter)
    {
        if (!String::valid(*iter))
        {
            return false;
        }
    }
    return true;
}


template<class String>
inline bool Foam::string::stripInvalid(std::string& str)
{
    if (!valid<String>(str))
    {
        size_type nValid = 0;
        iterator iter2 = str.begin();

        for
        (
            const_iterator iter1 = iter2;
            iter1 != const_cast<const std::string&>(str).end();
            iter1++
        )
        {
            const char c = *iter1;

            if (String::valid(c))
            {
                *iter2 = c;
                ++iter2;
                ++nValid;
            }
        }

        str.resize(nValid);

        return true;
    }

    return false;
}


template<class String>
inline bool Foam::string::meta(const std::string& str, const char quote)
{
    int escaped = 0;
    for (const_iterator iter = str.begin(); iter != str.end(); ++iter)
    {
        const char c = *iter;
        if (quote && c == quote)
        {
            escaped ^= 1;  // toggle state
        }
        else if (escaped)
        {
            escaped = 0;
        }
        else if (String::meta(c))
        {
            return true;
        }
    }
    return false;
}


template<class String>
inline Foam::string
Foam::string::quotemeta(const std::string& str, const char quote)
{
    if (!quote)
    {
        return str;
    }

    string sQuoted;
    sQuoted.reserve(2*str.length());

    int escaped = 0;
    for (const_iterator iter = str.begin(); iter != str.end(); ++iter)
    {
        const char c = *iter;
        if (c == quote)
        {
            escaped ^= 1;  // toggle state
        }
        else if (escaped)
        {
            escaped = 0;
        }
        else if (String::meta(c))
        {
            sQuoted += quote;
        }

        sQuoted += c;
    }

    sQuoted.resize(sQuoted.length());

    return sQuoted;
}


template<class String>
inline String Foam::string::validate(const std::string& str)
{
    string ss = str;
    stripInvalid<String>(ss);
    return ss;
}


inline bool Foam::string::match(const std::string& text) const
{
    return !compare(text);  // Always compare as literal string
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::string::operator()(const std::string& text) const
{
    return !compare(text);  // Always compare as literal string
}


inline Foam::string Foam::string::operator()
(
    const size_type i,
    const size_type n
) const
{
    return substr(i, n);
}


inline Foam::string Foam::string::operator()(const size_type n) const
{
    return substr(0, n);
}


inline unsigned Foam::string::hash::operator()
(
    const string& str,
    unsigned seed
) const
{
    return Hasher(str.data(), str.size(), seed);
}

// ************************************************************************* //

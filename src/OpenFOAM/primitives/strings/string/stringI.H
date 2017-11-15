/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline std::string::size_type Foam::string::find_ext() const
{
    const size_type i = find_last_of("./");

    if (i == npos || i == 0 || operator[](i) == '/')
    {
        return npos;
    }
    else
    {
        return i;
    }
}


inline bool Foam::string::hasExt() const
{
    return (find_ext() != npos);
}


inline bool Foam::string::removeExt()
{
    const size_type i = find_ext();

    if (i == npos)
    {
        return false;
    }
    else
    {
        this->resize(i);
        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::string::string()
{}


inline Foam::string::string(const std::string& str)
:
    std::string(str)
{}


inline Foam::string::string(const char* str)
:
    std::string(str)
{}


inline Foam::string::string(const char* str, const size_type len)
:
    std::string(str, len)
{}


inline Foam::string::string(const char c)
:
    std::string(1, c)
{}


inline Foam::string::string(const size_type len, const char c)
:
    std::string(len, c)
{}


inline Foam::string::string(std::string&& str)
:
    std::string(std::move(str))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class String>
inline bool Foam::string::valid(const std::string& str)
{
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
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
        size_type nChar = 0;
        iterator outIter = str.begin();

        for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
        {
            const char c = *iter;

            if (String::valid(c))
            {
                *outIter = c;
                ++outIter;
                ++nChar;
            }
        }

        str.resize(nChar);

        return true;
    }

    return false;
}


template<class String>
inline bool Foam::string::meta(const std::string& str, const char quote)
{
    int escaped = 0;
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
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
    sQuoted.reserve(2*str.size());

    int escaped = 0;
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
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

    sQuoted.resize(sQuoted.size());

    return sQuoted;
}


template<class String>
inline String Foam::string::validate(const std::string& str)
{
    String out;
    out.resize(str.size());

    size_type count = 0;
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        const char c = *iter;
        if (String::valid(c))
        {
            out[count++] = c;
        }
    }

    out.resize(count);

    return out;
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


inline unsigned Foam::string::hash::operator()
(
    const string& str,
    unsigned seed
) const
{
    return Hasher(str.data(), str.size(), seed);
}


// ************************************************************************* //

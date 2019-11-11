/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

inline std::string::size_type Foam::string::find_ext(const std::string& str)
{
    const auto i = str.find_last_of("./");

    if (i == npos || i == 0 || str[i] == '/')
    {
        return npos;
    }

    return i;
}


inline std::string::size_type Foam::string::find_ext() const
{
    return find_ext(*this);
}


inline bool Foam::string::hasPath() const
{
    return (npos != find('/'));
}


inline bool Foam::string::hasExt() const
{
    return (npos != find_ext());
}


inline bool Foam::string::removePath()
{
    const auto i = rfind('/');

    if (npos == i)
    {
        return false;
    }

    erase(0, i+1);
    return true;
}


inline bool Foam::string::removeExt()
{
    const auto i = find_ext();

    if (npos == i)
    {
        return false;
    }

    erase(i);
    return true;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::string::string(const std::string& str)
:
    std::string(str)
{}


inline Foam::string::string(std::string&& str)
:
    std::string(std::move(str))
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

        str.erase(nChar);

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

    sQuoted.shrink_to_fit();

    return sQuoted;
}


template<class String>
inline String Foam::string::validate(const std::string& str)
{
    String out;
    out.resize(str.size());

    size_type len = 0;
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        const char c = *iter;
        if (String::valid(c))
        {
            out[len] = c;
            ++len;
        }
    }

    out.erase(len);

    return out;
}


inline bool Foam::string::match(const std::string& text) const
{
    return !compare(text);  // Always compare as literal string
}


inline void Foam::string::swap(std::string& str)
{
    if (this != &str)
    {
        // Self-swap is a no-op
        std::string::swap(str);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::string::operator()(const std::string& text) const
{
    return !compare(text);  // Always compare as literal string
}


inline unsigned Foam::string::hash::operator()
(
    const std::string& str,
    unsigned seed
) const
{
    return Hasher(str.data(), str.size(), seed);
}


// ************************************************************************* //

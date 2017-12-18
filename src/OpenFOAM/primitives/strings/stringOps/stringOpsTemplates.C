/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include <cstdio>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NOTE: with C++11 could consider variadic templates for a more general
// sprintf implementation

template<class PrimitiveType>
Foam::word Foam::stringOps::name
(
    const char* fmt,
    const PrimitiveType& val
)
{
    word output;

    // snprintf with zero to find size (without '\0') required
    int n = ::snprintf(nullptr, 0, fmt, val);
    if (n > 0)
    {
        output.resize(n+1);
        char* buf = &(output[0]);

        // Print directly into buffer, no stripping desired
        n = ::snprintf(buf, n+1, fmt, val);
        output.resize(n);
    }

    return output;
}


template<class PrimitiveType>
Foam::word Foam::stringOps::name
(
    const std::string& fmt,
    const PrimitiveType& val
)
{
    return stringOps::name(fmt.c_str(), val);
}


template<class StringType>
Foam::SubStrings<StringType> Foam::stringOps::split
(
    const StringType& str,
    const char delim,
    const bool keepEmpty
)
{
    Foam::SubStrings<StringType> lst;
    if (str.empty() || !delim)
    {
        return lst;
    }

    lst.reserve(20);

    std::string::size_type beg = 0, end = 0;
    while ((end = str.find(delim, beg)) != std::string::npos)
    {
        if (keepEmpty || (beg < end))
        {
            lst.append(str.cbegin() + beg, str.cbegin() + end);
        }
        beg = end + 1;
    }

    // Trailing element
    if (keepEmpty ? (beg == str.size()) : (beg < str.size()))
    {
        lst.append(str.cbegin() + beg, str.cend());
    }

    return lst;
}


template<class StringType>
Foam::SubStrings<StringType> Foam::stringOps::split
(
    const StringType& str,
    const std::string& delim,
    const bool keepEmpty
)
{
    Foam::SubStrings<StringType> lst;
    if (str.empty() || delim.empty())
    {
        return lst;
    }

    lst.reserve(20);

    std::string::size_type beg = 0, end = 0;
    while ((end = str.find(delim, beg)) != std::string::npos)
    {
        if (keepEmpty || (beg < end))
        {
            lst.append(str.cbegin() + beg, str.cbegin() + end);
        }
        beg = end + delim.size();
    }

    // Trailing element
    if (keepEmpty ? (beg == str.size()) : (beg < str.size()))
    {
        lst.append(str.cbegin() + beg, str.cend());
    }

    return lst;
}


template<class StringType>
Foam::SubStrings<StringType> Foam::stringOps::splitAny
(
    const StringType& str,
    const std::string& delim
)
{
    Foam::SubStrings<StringType> lst;
    if (str.empty() || delim.empty())
    {
        return lst;
    }

    lst.reserve(20);

    for
    (
        std::string::size_type pos = 0;
        (pos = str.find_first_not_of(delim, pos)) != std::string::npos;
        /*nil*/
    )
    {
        const auto end = str.find_first_of(delim, pos);

        if (end == std::string::npos)
        {
            // Trailing element
            lst.append(str.cbegin() + pos, str.cend());
            break;
        }

        // Intermediate element
        lst.append(str.cbegin() + pos, str.cbegin() + end);

        pos = end + 1;
    }

    return lst;
}


template<class StringType>
Foam::SubStrings<StringType> Foam::stringOps::splitFixed
(
    const StringType& str,
    const std::string::size_type width,
    const std::string::size_type start
)
{
    Foam::SubStrings<StringType> lst;
    if (str.empty() || !width)
    {
        return lst;
    }

    const auto len = str.size();
    lst.reserve(1 + (len / width));

    for (std::string::size_type pos = start; pos < len; pos += width)
    {
        const auto end = (pos + width);

        if (end >= len)
        {
            // Trailing element
            lst.append(str.cbegin() + pos, str.cend());
            break;
        }

        lst.append(str.cbegin() + pos, str.cbegin() + end);
    }

    return lst;
}


template<class StringType>
Foam::SubStrings<StringType> Foam::stringOps::splitSpace
(
    const StringType& str
)
{
    return splitAny(str, "\t\n\v\f\r ");
}


// ************************************************************************* //

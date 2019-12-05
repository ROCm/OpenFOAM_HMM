/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline char* Foam::CStringList::stringCopy(char *dest, const std::string& str)
{
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        *dest = *iter;
        ++dest;
    }
    *(dest++) = '\0';

    return dest;
}


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline int Foam::CStringList::count(const char * const argv[])
{
    int n = 0;
    if (argv)
    {
        while (argv[n])
        {
            ++n;
        }
    }
    return n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::CStringList::CStringList()
:
    argc_(0),
    len_(0),
    argv_(nullptr),
    data_(nullptr)
{}


template<class StringType>
inline Foam::CStringList::CStringList(const UList<StringType>& input)
:
    CStringList()
{
    reset(input);
}


template<class StringType>
inline Foam::CStringList::CStringList(const SubStrings<StringType>& input)
:
    CStringList()
{
    reset(input);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::CStringList::~CStringList()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::CStringList::clear()
{
    argc_ = 0;
    len_  = 0;

    if (data_)
    {
        delete[] data_;
        data_ = nullptr;
    }
    if (argv_)
    {
        delete[] argv_;
        argv_ = nullptr;
    }
}


inline bool Foam::CStringList::empty() const noexcept
{
    return !argc_;
}


inline int Foam::CStringList::size() const noexcept
{
    return argc_;
}


inline size_t Foam::CStringList::length() const
{
    return len_;
}


inline char** Foam::CStringList::strings() const
{
    return argv_;
}


inline char** Foam::CStringList::strings(int start) const
{
    return &(argv_[start]);
}


inline const char* Foam::CStringList::data() const
{
    return data_;
}


template<class StringType>
inline int Foam::CStringList::reset(const UList<StringType>& input)
{
    return resetContent(input);
}


template<class StringType>
inline int Foam::CStringList::reset(const SubStrings<StringType>& input)
{
    return resetContent(input);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const char* Foam::CStringList::operator[](int i) const
{
    return argv_[i];
}


// ************************************************************************* //

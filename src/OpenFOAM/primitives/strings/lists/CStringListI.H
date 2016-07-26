/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

inline int Foam::CStringList::count(const char * const argv[])
{
    int nElem = 0;
    if (argv)
    {
        while (argv[nElem])
        {
            ++nElem;
        }
    }

    return nElem;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::CStringList::CStringList()
:
    argc_(0),
    len_(0),
    argv_(0),
    data_(0)
{}


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
        data_ = 0;
    }
    if (argv_)
    {
        delete[] argv_;
        argv_ = 0;
    }
}


inline int Foam::CStringList::size() const
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


inline char* Foam::CStringList::data() const
{
    return data_;
}


// ************************************************************************* //

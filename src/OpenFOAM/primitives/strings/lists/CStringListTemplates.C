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

#include "SubStrings.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ListType>
int Foam::CStringList::resetContent
(
    const ListType& input
)
{
    clear();

    if (input.empty())
    {
        // Special handling of an empty list
        argv_ = new char*[1];
        argv_[0] = nullptr;     // Final nullptr terminator
        return 0;
    }

    // Count overall required string length, including each trailing nul char
    for (const auto& str : input)
    {
        len_ += str.length() + 1;
    }
    --len_; // No final nul in overall count

    argv_ = new char*[input.size()+1];  // Extra +1 for terminating nullptr
    data_ = new char[len_+1];           // Extra +1 for terminating nul char

    argv_[0] = data_;   // Starts here

    for (const auto& str : input)
    {
        char *next = stringCopy(argv_[argc_], str);
        argv_[++argc_] = next;   // The start of next string
    }

    argv_[argc_] = nullptr;     // Final nullptr terminator

    return argc_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class StringType>
Foam::List<StringType>
Foam::CStringList::asList(int argc, const char * const argv[])
{
    List<StringType> lst(argc);

    for (int i=0; i < argc; ++i)
    {
        lst[i] = argv[i];
    }

    return lst;
}


template<class StringType>
Foam::List<StringType>
Foam::CStringList::asList(const char * const argv[])
{
    return asList<StringType>(count(argv), argv);
}


// ************************************************************************* //

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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class StringType>
Foam::CStringList::CStringList
(
    const UList<StringType>& input
)
:
    CStringList()
{
    reset(input);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class StringType>
void Foam::CStringList::reset
(
    const UList<StringType>& input
)
{
    clear();

    if (input.empty())
    {
        // Special handling of an empty list
        argv_ = new char*[1];
        argv_[0] = nullptr;     // Final nullptr terminator
        return;
    }

    // Count overall required string length, including each trailing nul char
    for (const auto& str : input)
    {
        len_ += str.size() + 1;
    }
    --len_; // No final nul in overall count

    argv_ = new char*[input.size()+1];  // Extra +1 for terminating nullptr
    data_ = new char[len_+1];           // Extra +1 for terminating nul char

    // Copy contents
    char* ptr = data_;
    for (const auto& str : input)
    {
        argv_[argc_++] = ptr;   // The start of this string

        for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
        {
            *(ptr++) = *iter;
        }
        *(ptr++) = '\0';
    }

    argv_[argc_] = nullptr;     // Final nullptr terminator
}


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

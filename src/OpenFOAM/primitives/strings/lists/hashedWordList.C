/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "hashedWordList.H"
#include "CStringList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hashedWordList::hashedWordList
(
    const label len,
    const char** array,
    bool unique
)
:
    wordList(len)
{
    for (label i=0; i < len; ++i)
    {
        wordList::operator[](i) = array[i];
    }

    rehash(unique);
}


Foam::hashedWordList::hashedWordList(const char** array, bool unique)
:
    hashedWordList(CStringList::count(array), array, unique)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hashedWordList::rehash() const
{
    lookup_.clear();

    const wordUList& list = *this;
    const label len = list.size();

    for (label i=0; i < len; ++i)
    {
        lookup_.insert(list[i], i);
    }
}


void Foam::hashedWordList::uniq()
{
    lookup_.clear();

    wordList& list = *this;
    const label len = list.size();

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        word& item = list[i];

        if (lookup_.insert(item, i))
        {
            if (count != i)
            {
                list[count] = std::move(item);
            }
            ++count;
        }
    }

    list.resize(count);
}


// ************************************************************************* //

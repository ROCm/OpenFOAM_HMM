/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hashedWordList::hashedWordList
(
    const label count,
    const char** lst,
    const bool removeDuplicates
)
:
    List<word>(count)
{
    forAll(*this, i)
    {
        List<word>::operator[](i) = lst[i];
    }

    rehash(removeDuplicates);
}


Foam::hashedWordList::hashedWordList
(
    const char** lst,
    const bool removeDuplicates
)
{
    // Determine the number of entries
    label count = 0;
    for (unsigned i = 0; lst[i] && *(lst[i]); ++i)
    {
        ++count;
    }

    List<word>::setSize(count);
    forAll(*this, i)
    {
        List<word>::operator[](i) = lst[i];
    }

    rehash(removeDuplicates);
}


Foam::hashedWordList::hashedWordList(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::hashedWordList::transfer
(
    List<word>& lst,
    const bool removeDuplicates
)
{
    List<word>::transfer(lst);
    rehash(removeDuplicates);
}


void Foam::hashedWordList::rehash() const
{
    indices_.clear();

    forAll(*this, i)
    {
        indices_.insert(List<word>::operator[](i), i);
    }
}


void Foam::hashedWordList::uniq()
{
    indices_.clear();

    label nElem = 0;
    forAll(*this, i)
    {
        const word& item = List<word>::operator[](i);

        if (indices_.insert(item, nElem))
        {
            if (nElem != i)
            {
                List<word>::operator[](nElem) = item;
            }
            ++nElem;
        }
    }

    List<word>::setSize(nElem);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, hashedWordList& lst)
{
    is  >> static_cast<List<word>&>(lst);
    lst.rehash();

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const hashedWordList& lst)
{
    os  << static_cast<const UList<word>&>(lst);
    return os;
}


// ************************************************************************* //

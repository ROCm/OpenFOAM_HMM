/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "UPtrList.H"
#include "PtrList.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::UPtrList<T>::UPtrList(PtrList<T>& list)
:
    ptrs_(list.ptrs_) // shallow copy (via const reference)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::UPtrList<T>::reorder(const labelUList& oldToNew)
{
    const label len = this->size();

    if (oldToNew.size() != len)
    {
        FatalErrorInFunction
            << "Size of map (" << oldToNew.size()
            << ") not equal to list size (" << len
            << ") for type " << typeid(T).name() << nl
            << abort(FatalError);
    }

    Detail::PtrListDetail<T> newList(len);

    for (label i=0; i<len; ++i)
    {
        const label idx = oldToNew[i];

        if (idx < 0 || idx >= len)
        {
            FatalErrorInFunction
                << "Illegal index " << idx << nl
                << "Valid indices are [0," << len << ") for type "
                << typeid(T).name() << nl
                << abort(FatalError);
        }

        if (newList[idx])
        {
            FatalErrorInFunction
                << "reorder map is not unique; element " << idx
                << " already set for type " << typeid(T).name()
                << abort(FatalError);
        }
        newList[idx] = ptrs_[i];
    }

    // Verify that all pointers were indeed set
    for (label i=0; i<len; ++i)
    {
        if (!newList[i])
        {
            FatalErrorInFunction
                << "Element " << i << " not set after reordering." << nl
                << abort(FatalError);
        }
    }

    ptrs_.transfer(newList);
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const UPtrList<T>& list)
{
    return list.ptrs_.write(os);
}


// ************************************************************************* //

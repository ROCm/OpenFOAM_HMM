/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "PackedBoolList.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PackedBoolList::PackedBoolList(Istream& is)
:
    PackedList<1>()
{
    is  >> *this;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList> Foam::PackedBoolList::used() const
{
    labelList lst(this->count());

    if (lst.size())
    {
        label nElem = 0;

        forAll(*this, elemI)
        {
            if (get(elemI))
            {
                lst[nElem++] = elemI;
            }
        }

        lst.setSize(nElem);
    }

    return lst.xfer();
}


template<class LabelListType>
Foam::label Foam::PackedBoolList::setFromIndices(const LabelListType& indices)
{
    // no better information, just guess something about the size
    reserve(indices.size());

    label cnt = 0;
    forAll(indices, elemI)
    {
        if (set(indices[elemI]))
        {
            ++cnt;
        }
    }

    return cnt;
}


template<class LabelListType>
Foam::label Foam::PackedBoolList::unsetFromIndices(const LabelListType& indices)
{
    label cnt = 0;
    forAll(indices, elemI)
    {
        if (unset(indices[elemI]))
        {
            ++cnt;
        }
    }

    return cnt;
}


Foam::label Foam::PackedBoolList::set(const UList<label>& indices)
{
    return setFromIndices(indices);
}


Foam::label Foam::PackedBoolList::set(const UIndirectList<label>& indices)
{
    return setFromIndices(indices);
}


Foam::label Foam::PackedBoolList::unset(const UList<label>& indices)
{
    return unsetFromIndices(indices);
}


Foam::label Foam::PackedBoolList::unset(const UIndirectList<label>& indices)
{
    return unsetFromIndices(indices);
}


void Foam::PackedBoolList::modulo(const PackedList<1>& lst)
{
    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    const label len = min(lhs.size(), rhs.size());

    for (label i=0; i < len; ++i)
    {
        lhs[i] &= ~rhs[i];
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


Foam::PackedBoolList&
Foam::PackedBoolList::operator=(const UList<bool>& lst)
{
    this->setSize(lst.size());

    forAll(*this, elemI)
    {
        set(elemI, lst[elemI]);
    }

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator=(const UList<label>& indices)
{
    clear();
    set(indices);

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator=(const UIndirectList<label>& indices)
{
    clear();
    set(indices);

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator|=(const PackedList<1>& lst)
{
    // extend addressable area if needed
    if (this->size() < lst.size())
    {
        this->resize(lst.size());
    }

    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    forAll(rhs, i)
    {
        lhs[i] |= rhs[i];
    }

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator&=(const PackedList<1>& lst)
{
    // shrink addressable area if needed
    if (this->size() > lst.size())
    {
        this->resize(lst.size());
    }

    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    forAll(lhs, i)
    {
        lhs[i] &= rhs[i];
    }

    // trim to bits actually used
    this->trim();

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator^=(const PackedList<1>& lst)
{
    // extend addressable area if needed
    if (this->size() < lst.size())
    {
        this->resize(lst.size());
    }

    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    forAll(rhs, i)
    {
        lhs[i] ^= rhs[i];
    }

    return *this;
}


// * * * * * * * * * * * * * *  Global Operators * * * * * * * * * * * * * * //

Foam::PackedBoolList Foam::operator|
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
)
{
    PackedBoolList result(lst1);
    result |= lst2;
    return result;
}


Foam::PackedBoolList Foam::operator&
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
)
{
    PackedBoolList result(lst1);
    result &= lst2;
    return result;
}


Foam::PackedBoolList Foam::operator^
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
)
{
    PackedBoolList result(lst1);
    result ^= lst2;
    return result;
}


// ************************************************************************* //

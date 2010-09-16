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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::PackedBoolList::bitorPrepare
(
    const PackedList<1>& lst,
    label& maxPackLen
)
{
    const StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    const label packLen1 = this->packedLength();
    const label packLen2 = lst.packedLength();


    // check how the lists interact and if bit trimming is needed
    bool needTrim = false;
    maxPackLen = packLen1;

    if (packLen1 == packLen2)
    {
        // identical packed lengths - only resize if absolutely necessary
        if
        (
            this->size() != lst.size()
         && maxPackLen
         && rhs[maxPackLen-1] > lhs[maxPackLen-1]
        )
        {
            // second list has a higher bit set
            // extend addressable area and use trim
            resize(lst.size());
            needTrim = true;
        }
    }
    else if (packLen2 < packLen1)
    {
        // second list is shorter, this limits the or
        maxPackLen = packLen2;
    }
    else
    {
        // second list is longer, find the highest bit set
        for (label storeI = packLen1; storeI < packLen2; ++storeI)
        {
            if (rhs[storeI])
            {
                maxPackLen = storeI+1;
            }
        }

        // the upper limit moved - resize for full coverage and trim later
        if (maxPackLen > packLen1)
        {
            resize(maxPackLen * packing());
            needTrim = true;
        }
    }

    return needTrim;
}



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
Foam::label Foam::PackedBoolList::setIndices(const LabelListType& indices)
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
Foam::label Foam::PackedBoolList::unsetIndices(const LabelListType& indices)
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
    return setIndices(indices);
}


Foam::label Foam::PackedBoolList::set(const UIndirectList<label>& indices)
{
    return setIndices(indices);
}


Foam::label Foam::PackedBoolList::unset(const UList<label>& indices)
{
    return unsetIndices(indices);
}


Foam::label Foam::PackedBoolList::unset(const UIndirectList<label>& indices)
{
    return unsetIndices(indices);
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

    const label len = this->packedLength();

    for (label i=0; i < len; ++i)
    {
        lhs[i] &= rhs[i];
    }

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator^=(const PackedList<1>& lst)
{
    // extend addressable area if needed, return maximum size possible
    label len = 0;
    const bool needTrim = bitorPrepare(lst, len);

    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    for (label i=0; i < len; ++i)
    {
        lhs[i] ^= rhs[i];
    }

    if (needTrim)
    {
        trim();
    }

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator|=(const PackedList<1>& lst)
{
    // extend addressable area if needed, return maximum size possible
    label len = 0;
    const bool needTrim = bitorPrepare(lst, len);

    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    for (label i=0; i < len; ++i)
    {
        lhs[i] |= rhs[i];
    }

    if (needTrim)
    {
        trim();
    }

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator-=(const PackedList<1>& lst)
{
    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    // overlapping storage size
    const label len = min(this->packedLength(), lst.packedLength());

    for (label i=0; i < len; ++i)
    {
        lhs[i] &= ~rhs[i];
    }

    return *this;
}


// * * * * * * * * * * * * * *  Global Operators * * * * * * * * * * * * * * //

Foam::PackedBoolList Foam::operator&
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
)
{
    PackedBoolList result(lst1);
    result &= lst2;

    // trim to bits actually used
    result.trim();

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

    // trim to bits actually used
    result.trim();

    return result;
}


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


// ************************************************************************* //

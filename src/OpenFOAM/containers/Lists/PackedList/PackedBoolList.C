/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    const PackedBoolList& lst,
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


template<class LabelListType>
void Foam::PackedBoolList::setIndices(const LabelListType& indices)
{
    const label len = indices.size();

    // No better information, just guess something from the size
    reserve(len);

    for (label i = 0; i < len; ++i)
    {
        set(indices[i]);
    }
}


template<class LabelListType>
void Foam::PackedBoolList::unsetIndices(const LabelListType& indices)
{
    const label len = indices.size();
    for (label i = 0; i < len; ++i)
    {
        unset(indices[i]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PackedBoolList::PackedBoolList(Istream& is)
:
    PackedList<1>()
{
    is  >> *this;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PackedBoolList::set(const PackedBoolList& lst)
{
    // extend addressable area if needed, return maximum size possible
    label len = 0;
    const bool needTrim = bitorPrepare(lst, len);

    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    for (label i = 0; i < len; ++i)
    {
        lhs[i] |= rhs[i];
    }

    if (needTrim)
    {
        trim();
    }
}


void Foam::PackedBoolList::unset(const PackedBoolList& lst)
{
    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    // overlapping storage size
    const label len = min(this->packedLength(), lst.packedLength());

    for (label i = 0; i < len; ++i)
    {
        lhs[i] &= ~rhs[i];
    }
}


void Foam::PackedBoolList::setMany(const labelUList& indices)
{
    setIndices(indices);
}


void Foam::PackedBoolList::setMany(const labelUIndList& indices)
{
    setIndices(indices);
}


void Foam::PackedBoolList::unsetMany(const labelUList& indices)
{
    unsetIndices(indices);
}


void Foam::PackedBoolList::unsetMany(const labelUIndList& indices)
{
    unsetIndices(indices);
}


Foam::labelList Foam::PackedBoolList::used() const
{
    // Number of used (set) entries
    const label cnt = this->count();

    labelList lst(cnt);

    if (cnt)
    {
        // The length of the input list
        const label len = this->size();

        for (label i=0, usedi=0; (i < len && usedi < cnt); ++i)
        {
            if (test(i))
            {
                lst[usedi++] = i;
            }
        }
    }

    return lst;
}


// ************************************************************************* //

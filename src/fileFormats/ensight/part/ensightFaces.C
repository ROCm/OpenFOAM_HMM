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

#include "ensightFaces.H"
#include "error.H"
#include "polyMesh.H"
#include "ListOps.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::ensightFaces::nTypes = 3;

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::ensightFaces::elemType,
        3
    >::names[] = { "tria3", "quad4", "nsided" };
}

const Foam::NamedEnum<Foam::ensightFaces::elemType, 3>
    Foam::ensightFaces::elemEnum;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// only used in this file-scope
inline Foam::ensightFaces::elemType
Foam::ensightFaces::whatType(const face& f)
{
    return
    (
        f.size() == 3
      ? Foam::ensightFaces::elemType::TRIA3
      : f.size() == 4
      ? Foam::ensightFaces::elemType::QUAD4
      : Foam::ensightFaces::elemType::NSIDED
    );
}


// only used in this file-scope
inline void Foam::ensightFaces::add
(
    const face& f,
    const label id,
    const bool flip
)
{
    const enum elemType what = whatType(f);

    // linear addressing:
    const label index = offset(what) + sizes_[what]++;

    address_[index] = id;
    if (flipMap_.size())
    {
        flipMap_[index] = flip;
    }
}


void Foam::ensightFaces::resizeAll()
{
    // overall required size
    label n = 0;
    forAll(sizes_, typei)
    {
        n += sizes_[typei];
    }
    address_.setSize(n, Zero);

    // assign corresponding sub-lists
    n = 0;
    forAll(sizes_, typei)
    {
        deleteDemandDrivenData(lists_[typei]);

        lists_[typei] = new SubList<label>(address_, sizes_[typei], n);

        n += sizes_[typei];
    }

    // normally assume no flipMap
    flipMap_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightFaces::ensightFaces(label partIndex)
:
    index_(partIndex),
    address_(),
    flipMap_(),
    sizes_(Zero),
    lists_()
{
    // Ensure sub-lists are properly initialized to nullptr
    forAll(lists_, typei)
    {
        lists_[typei] = nullptr;
    }

    resizeAll(); // adjust allocation
}


Foam::ensightFaces::ensightFaces(const ensightFaces& obj)
:
    index_(obj.index_),
    address_(obj.address_),
    flipMap_(obj.flipMap_),
    sizes_(),
    lists_()
{
    // Ensure sub-lists are properly initialized to nullptr
    forAll(lists_, typei)
    {
        lists_[typei] = nullptr;
    }

    // Total (reduced) sizes
    FixedList<label, 3> totSizes = obj.sizes_;

    // Local sizes
    this->sizes_ = obj.sizes();

    resizeAll(); // adjust allocation

    // Restore total (reduced) sizes
    this->sizes_ = totSizes;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightFaces::~ensightFaces()
{
    forAll(lists_, typei)
    {
        deleteDemandDrivenData(lists_[typei]);
    }
    address_.clear();
    flipMap_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::label, 3> Foam::ensightFaces::sizes() const
{
    FixedList<label, 3> count;
    forAll(lists_, typei)
    {
        count[typei] = lists_[typei]->size();
    }

    return count;
}


Foam::label Foam::ensightFaces::offset(const enum elemType what) const
{
    label n = 0;
    for (label typei = 0; typei < label(what); ++typei)
    {
        n += lists_[typei]->size();
    }

    return n;
}


Foam::label Foam::ensightFaces::total() const
{
    label n = 0;
    forAll(sizes_, typei)
    {
        n += sizes_[typei];
    }
    return n;
}


void Foam::ensightFaces::clear()
{
    sizes_ = Zero;  // reset sizes
    resizeAll();
}


void Foam::ensightFaces::reduce()
{
    forAll(sizes_, typei)
    {
        sizes_[typei] = lists_[typei]->size();
        Foam::reduce(sizes_[typei], sumOp<label>());
    }
}


void Foam::ensightFaces::sort()
{
    if (flipMap_.size() == address_.size())
    {
        // sort flip map too

        labelList order;
        label start = 0;

        forAll(lists_, typei)
        {
            SubList<label>& idLst = *(lists_[typei]);
            const label sz = idLst.size();

            if (sz)
            {
                SubList<bool> flip(flipMap_, sz, start);
                start += sz; // for next sub-list

                Foam::sortedOrder(idLst, order);

                idLst = reorder<labelList>(order, idLst);
                flip  = reorder<boolList>(order,  flip);
            }
        }
    }
    else
    {
        // no flip-maps, simpler to sort
        forAll(lists_, typei)
        {
            Foam::sort(*(lists_[typei]));
        }
        flipMap_.clear();  // for safety
    }
}


void Foam::ensightFaces::classify(const faceList& faces)
{
    const label sz = faces.size();

    // Count the shapes
    // Can avoid double looping, but only at the expense of allocation

    sizes_ = Zero;  // reset sizes
    for (label listi = 0; listi < sz; ++listi)
    {
        const enum elemType what = whatType(faces[listi]);
        sizes_[what]++;
    }

    resizeAll();    // adjust allocation
    sizes_ = Zero;  // reset sizes

    // Assign face-id per shape type
    for (label listi = 0; listi < sz; ++listi)
    {
        add(faces[listi], listi);
    }
}


void Foam::ensightFaces::classify
(
    const faceList& faces,
    const labelUList& addressing,
    const boolList& flipMap,
    const PackedBoolList& exclude
)
{
    // Note: Since PackedList::operator[] returns zero (false) for out-of-range
    // indices, can skip our own bounds checking here.

    const label sz = addressing.size();
    const bool useFlip = (addressing.size() == flipMap.size());

    // Count the shapes
    // Can avoid double looping, but only at the expense of allocation

    sizes_ = Zero;  // reset sizes
    for (label listi = 0; listi < sz; ++listi)
    {
        const label faceId = addressing[listi];

        if (!exclude[faceId])
        {
            const enum elemType what = whatType(faces[faceId]);
            sizes_[what]++;
        }
    }

    resizeAll();    // adjust allocation
    sizes_ = Zero;  // reset sizes

    if (useFlip)
    {
        flipMap_.setSize(address_.size(), false);
        flipMap_ = false;
    }

    // Assign face-id per shape type
    for (label listi = 0; listi < sz; ++listi)
    {
        const label faceId = addressing[listi];
        const bool flip = useFlip && flipMap[listi];

        if (!exclude[faceId])
        {
            add(faces[faceId], faceId, flip);
        }
    }
}

// ************************************************************************* //

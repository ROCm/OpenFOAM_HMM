/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::ensightFaces::nTypes = 3;

const char* Foam::ensightFaces::elemNames[3] =
    { "tria3", "quad4", "nsided" };


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
        slices_[typei].setStart(n);
        slices_[typei].setSize(sizes_[typei]);

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
    slices_(),
    sizes_(Zero)
{
    resizeAll(); // adjust allocation/sizing
}


Foam::ensightFaces::ensightFaces(const ensightFaces& obj)
:
    index_(obj.index_),
    address_(obj.address_),
    flipMap_(obj.flipMap_),
    slices_(),
    sizes_()
{
    // Save the total (reduced) sizes
    FixedList<label, 3> totSizes = obj.sizes_;

    // Need local sizes for the resize operation
    this->sizes_ = obj.sizes();

    resizeAll(); // adjust allocation/sizing

    // Restore total (reduced) sizes
    this->sizes_ = totSizes;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightFaces::~ensightFaces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::label, 3> Foam::ensightFaces::sizes() const
{
    FixedList<label, 3> count;
    forAll(slices_, typei)
    {
        count[typei] = slices_[typei].size();
    }

    return count;
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
    // No listCombineGather, listCombineScatter for FixedList
    forAll(sizes_, typei)
    {
        sizes_[typei] = slices_[typei].size();
        Foam::reduce(sizes_[typei], sumOp<label>());
    }
}


void Foam::ensightFaces::sort()
{
    if (flipMap_.size() == address_.size())
    {
        // Must sort flip map as well
        labelList order;

        forAll(slices_, typei)
        {
            if (slices_[typei].size())
            {
                SubList<label> idLst(address_, slices_[typei]);
                SubList<bool>  flip(flipMap_, slices_[typei]);

                Foam::sortedOrder(idLst, order);

                idLst = reorder<labelList>(order, idLst);
                flip  = reorder<boolList>(order,  flip);
            }
        }
    }
    else
    {
        // no flip-maps, simpler to sort
        forAll(slices_, typei)
        {
            if (slices_[typei].size())
            {
                SubList<label> idLst(address_, slices_[typei]);
                Foam::sort(idLst);
            }
        }

        flipMap_.clear();  // for extra safety
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
    sizes_ = Zero;  // reset sizes - use for local indexing here

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
    sizes_ = Zero;  // reset sizes - use for local indexing here

    if (useFlip)
    {
        flipMap_.setSize(address_.size(), false);
        flipMap_ = false;
    }

    // Assign face-id per shape type
    for (label listi = 0; listi < sz; ++listi)
    {
        const label faceId = addressing[listi];
        const bool  doFlip = useFlip && flipMap[listi];

        if (!exclude[faceId])
        {
            add(faces[faceId], faceId, doFlip);
        }
    }
}

// ************************************************************************* //

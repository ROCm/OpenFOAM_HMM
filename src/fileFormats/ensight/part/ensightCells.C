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

#include "ensightCells.H"
#include "error.H"
#include "polyMesh.H"
#include "cellModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::ensightCells::nTypes = 5;

const char* Foam::ensightCells::elemNames[5] =
    { "tetra4", "pyramid5", "penta6", "hexa8", "nfaced" };


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightCells::resizeAll()
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightCells::ensightCells(const label partIndex)
:
    index_(partIndex),
    address_(),
    slices_(),
    sizes_(Zero)
{
    resizeAll(); // adjust allocation/sizing
}


Foam::ensightCells::ensightCells(const ensightCells& obj)
:
    index_(obj.index_),
    address_(obj.address_),
    slices_(),
    sizes_()
{
    // Save the total (reduced) sizes
    FixedList<label, 5> totSizes = obj.sizes_;

    // Need local sizes for the resize operation
    this->sizes_ = obj.sizes();

    resizeAll(); // adjust allocation/sizing

    // Restore total (reduced) sizes
    this->sizes_ = totSizes;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightCells::~ensightCells()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::label, 5> Foam::ensightCells::sizes() const
{
    FixedList<label, 5> count;
    forAll(slices_, typei)
    {
        count[typei] = slices_[typei].size();
    }

    return count;
}


Foam::label Foam::ensightCells::total() const
{
    label n = 0;
    forAll(sizes_, typei)
    {
        n += sizes_[typei];
    }
    return n;
}


void Foam::ensightCells::clear()
{
    sizes_ = Zero;  // reset sizes
    resizeAll();
}


void Foam::ensightCells::reduce()
{
    // No listCombineGather, listCombineScatter for FixedList
    forAll(sizes_, typei)
    {
        sizes_[typei] = slices_[typei].size();
        Foam::reduce(sizes_[typei], sumOp<label>());
    }
}


void Foam::ensightCells::sort()
{
    forAll(slices_, typei)
    {
        if (slices_[typei].size())
        {
            SubList<label> idLst(address_, slices_[typei]);
            Foam::sort(idLst);
        }
    }
}


void Foam::ensightCells::classify
(
    const polyMesh& mesh,
    const labelUList& addressing
)
{
    // References to cell shape models
    const cellModel& tet   = cellModel::ref(cellModel::TET);
    const cellModel& pyr   = cellModel::ref(cellModel::PYR);
    const cellModel& prism = cellModel::ref(cellModel::PRISM);
    const cellModel& hex   = cellModel::ref(cellModel::HEX);

    const cellShapeList& shapes = mesh.cellShapes();

    const bool indirect = notNull(addressing);
    const label sz = indirect ? addressing.size() : mesh.nCells();

    // Count the shapes
    // Can avoid double looping, but only at the expense of allocation

    sizes_ = Zero;  // reset sizes
    for (label listi = 0; listi < sz; ++listi)
    {
        const label id = indirect ? addressing[listi] : listi;
        const cellModel& model = shapes[id].model();

        enum elemType what = NFACED;
        if (model == tet)
        {
            what = TETRA4;
        }
        else if (model == pyr)
        {
            what = PYRAMID5;
        }
        else if (model == prism)
        {
            what = PENTA6;
        }
        else if (model == hex)
        {
            what = HEXA8;
        }

        sizes_[what]++;
    }

    resizeAll();    // adjust allocation
    sizes_ = Zero;  // reset sizes - use for local indexing here

    // Assign cell-id per shape type
    for (label listi = 0; listi < sz; ++listi)
    {
        const label id = indirect ? addressing[listi] : listi;
        const cellModel& model = shapes[id].model();

        enum elemType what = NFACED;
        if (model == tet)
        {
            what = TETRA4;
        }
        else if (model == pyr)
        {
            what = PYRAMID5;
        }
        else if (model == prism)
        {
            what = PENTA6;
        }
        else if (model == hex)
        {
            what = HEXA8;
        }

        // eg, the processor local cellId
        labelUList slice = address_[slices_[what]];

        slice[sizes_[what]] = id;
        sizes_[what]++;
    }
}


// ************************************************************************* //

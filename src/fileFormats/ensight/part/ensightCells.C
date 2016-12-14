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

#include "ensightCells.H"
#include "error.H"
#include "polyMesh.H"
#include "cellModeller.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::ensightCells::nTypes = 5;

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::ensightCells::elemType,
        5
    >::names[] = { "tetra4", "pyramid5", "penta6", "hexa8", "nfaced" };
}

const Foam::NamedEnum<Foam::ensightCells::elemType, 5>
    Foam::ensightCells::elemEnum;


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
        deleteDemandDrivenData(lists_[typei]);

        lists_[typei] = new SubList<label>(address_, sizes_[typei], n);

        n += sizes_[typei];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightCells::ensightCells(const label partIndex)
:
    index_(partIndex),
    address_(),
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


Foam::ensightCells::ensightCells(const ensightCells& obj)
:
    index_(obj.index_),
    address_(obj.address_),
    sizes_(),
    lists_()
{
    // Ensure sub-lists are properly initialized to nullptr
    forAll(lists_, typei)
    {
        lists_[typei] = nullptr;
    }

    // Total (reduced) sizes
    FixedList<label, 5> totSizes = obj.sizes_;

    // Local sizes
    this->sizes_ = obj.sizes();

    resizeAll(); // adjust allocation

    // Restore total (reduced) sizes
    this->sizes_ = totSizes;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightCells::~ensightCells()
{
    forAll(lists_, typei)
    {
        deleteDemandDrivenData(lists_[typei]);
    }
    address_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::label, 5> Foam::ensightCells::sizes() const
{
    FixedList<label, 5> count;
    forAll(lists_, typei)
    {
        count[typei] = lists_[typei]->size();
    }

    return count;
}


Foam::label Foam::ensightCells::offset(const enum elemType what) const
{
    label n = 0;
    for (label typei = 0; typei < label(what); ++typei)
    {
        n += lists_[typei]->size();
    }

    return n;
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
    forAll(sizes_, typei)
    {
        sizes_[typei] = lists_[typei]->size();
        Foam::reduce(sizes_[typei], sumOp<label>());
    }
}


void Foam::ensightCells::sort()
{
    forAll(lists_, typei)
    {
        Foam::sort(*(lists_[typei]));
    }
}


void Foam::ensightCells::classify
(
    const polyMesh& mesh,
    const labelUList& addressing
)
{
    // References to cell shape models
    const cellModel& tet   = *(cellModeller::lookup("tet"));
    const cellModel& pyr   = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& hex   = *(cellModeller::lookup("hex"));

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
    sizes_ = Zero;  // reset sizes

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
        lists_[what]->operator[](sizes_[what]++) = id;
    }
}


// ************************************************************************* //

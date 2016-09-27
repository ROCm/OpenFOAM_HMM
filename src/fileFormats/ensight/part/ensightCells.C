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

const Foam::NamedEnum<Foam::ensightCells::elemType,5>
    Foam::ensightCells::elemEnum;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::label Foam::ensightCells::offset
(
    const enum elemType what,
    const label i
) const
{
    label n = i;
    for (label typeI = 0; typeI < label(what); ++typeI)
    {
        n += sizes_[typeI];
    }

    return n;
}


void Foam::ensightCells::allocate()
{
    // overall required size
    label n = 0;
    forAll(sizes_, typeI)
    {
        n += sizes_[typeI];
    }
    address_.setSize(n, 0);

    // assign corresponding sub-lists
    n = 0;
    forAll(sizes_, typeI)
    {
        deleteDemandDrivenData(lists_[typeI]);

        lists_[typeI] = new SubList<label>(address_, sizes_[typeI], n);

        n += sizes_[typeI];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightCells::ensightCells(const label partIndex)
:
    index_(partIndex),
    address_(),
    sizes_(0),
    lists_()
{
    // ensure sub-lists are properly initialized to nullptr
    forAll(lists_, typeI)
    {
        lists_[typeI] = nullptr;
    }

    clear();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightCells::~ensightCells()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::label, 5> Foam::ensightCells::sizes() const
{
    FixedList<label, 5> count;
    forAll(lists_, typeI)
    {
        count[typeI] = lists_[typeI] ? lists_[typeI]->size() : 0;
    }

    return count;
}


Foam::label Foam::ensightCells::total() const
{
    label n = 0;
    forAll(sizes_, typeI)
    {
        n += sizes_[typeI];
    }
    return n;
}


void Foam::ensightCells::clear()
{
    sizes_ = 0;

    forAll(lists_, typeI)
    {
        deleteDemandDrivenData(lists_[typeI]);
    }
    address_.clear();
}


void Foam::ensightCells::reduce()
{
    forAll(sizes_, typeI)
    {
        sizes_[typeI] = lists_[typeI] ? lists_[typeI]->size() : 0;
        Foam::reduce(sizes_[typeI], sumOp<label>());
    }
}


void Foam::ensightCells::sort()
{
    forAll(lists_, typeI)
    {
        if (lists_[typeI])
        {
            Foam::sort(*(lists_[typeI]));
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
    const cellModel& tet   = *(cellModeller::lookup("tet"));
    const cellModel& pyr   = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& hex   = *(cellModeller::lookup("hex"));

    const cellShapeList& shapes = mesh.cellShapes();

    const bool indirect = notNull(addressing);
    const label sz = indirect ? addressing.size() : mesh.nCells();

    // Count the shapes
    // Can avoid double looping, but only at the expense of allocation

    sizes_ = 0; // reset sizes
    for (label listI = 0; listI < sz; ++listI)
    {
        const label id = indirect ? addressing[listI] : listI;
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

    allocate();
    sizes_ = 0;   // reset sizes

    // Assign cell-id per shape type
    for (label listI = 0; listI < sz; ++listI)
    {
        const label id = indirect ? addressing[listI] : listI;
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


Foam::label Foam::ensightCells::offset(const enum elemType what) const
{
    return offset(what, 0);
}


// ************************************************************************* //

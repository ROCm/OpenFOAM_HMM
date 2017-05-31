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

#include "polyMesh.H"
#include "foamVtkCells.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkCells::foamVtkCells
(
    const contentType output,
    const bool decompose
)
:
    foamVtuSizing(),
    output_(output),
    decomposeRequest_(decompose),
    cellTypes_(),
    vertLabels_(),
    vertOffset_(),
    faceLabels_(),
    faceOffset_(),
    maps_()
{}


Foam::foamVtkCells::foamVtkCells
(
    const polyMesh& mesh,
    const contentType output,
    const bool decompose
)
:
    foamVtkCells(output, decompose)
{
    reset(mesh);
}


Foam::foamVtkCells::foamVtkCells
(
    const foamVtkOutput::outputOptions outOpts,
    const bool decompose
)
:
    foamVtkCells
    (
        (outOpts.legacy() ? contentType::LEGACY : contentType::XML),
        decompose
    )
{}


Foam::foamVtkCells::foamVtkCells
(
    const polyMesh& mesh,
    const foamVtkOutput::outputOptions outOpts,
    const bool decompose
)
:
    foamVtkCells(outOpts, decompose)
{
    reset(mesh);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamVtkCells::~foamVtkCells()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::foamVtkCells::clear()
{
    foamVtuSizing::clear();
    cellTypes_.clear();
    vertLabels_.clear();
    vertOffset_.clear();
    faceLabels_.clear();
    faceOffset_.clear();

    maps_.clear();
}


void Foam::foamVtkCells::reset(const polyMesh& mesh)
{
    foamVtuSizing::reset(mesh, decomposeRequest_);

    cellTypes_.setSize(nFieldCells());
    vertLabels_.setSize(slotSize(output_, slotType::CELLS));
    vertOffset_.setSize(slotSize(output_, slotType::CELLS_OFFSETS));
    faceLabels_.setSize(slotSize(output_, slotType::FACES));
    faceOffset_.setSize(slotSize(output_, slotType::FACES_OFFSETS));

    switch (output_)
    {
        case contentType::LEGACY:
            populateLegacy
            (
                mesh,
                cellTypes_,
                vertLabels_,
                maps_
            );
            break;
        case contentType::XML:
            populateXml
            (
                mesh,
                cellTypes_,
                vertLabels_,
                vertOffset_,
                faceLabels_,
                faceOffset_,
                maps_
            );
            break;
        case contentType::INTERNAL:
            populateInternal
            (
                mesh,
                cellTypes_,
                vertLabels_,
                vertOffset_,
                faceLabels_,
                faceOffset_,
                maps_
            );
            break;
    }
}


void Foam::foamVtkCells::reset
(
    const polyMesh& mesh,
    const enum contentType output,
    const bool decompose
)
{
    output_ = output;
    decomposeRequest_ = decompose;

    reset(mesh);
}


void Foam::foamVtkCells::renumberCells(const UList<label>& mapping)
{
    maps_.renumberCells(mapping);
}


void Foam::foamVtkCells::renumberPoints(const UList<label>& mapping)
{
    maps_.renumberPoints(mapping);
}


// ************************************************************************* //

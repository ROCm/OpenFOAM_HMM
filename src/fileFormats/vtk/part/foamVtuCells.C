/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
#include "foamVtuCells.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::vtuCells::vtuCells
(
    const contentType output,
    const bool decompose
)
:
    vtk::vtuSizing(),
    output_(output),
    decomposeRequest_(decompose),
    cellTypes_(),
    vertLabels_(),
    vertOffset_(),
    faceLabels_(),
    faceOffset_(),
    maps_()
{}


Foam::vtk::vtuCells::vtuCells
(
    const polyMesh& mesh,
    const contentType output,
    const bool decompose
)
:
    vtuCells(output, decompose)
{
    reset(mesh);
}


Foam::vtk::vtuCells::vtuCells
(
    const vtk::outputOptions opts,
    const bool decompose
)
:
    vtuCells
    (
        (opts.legacy() ? contentType::LEGACY : contentType::XML),
        decompose
    )
{}


Foam::vtk::vtuCells::vtuCells
(
    const polyMesh& mesh,
    const vtk::outputOptions opts,
    const bool decompose
)
:
    vtuCells(opts, decompose)
{
    reset(mesh);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::vtuCells::clear()
{
    vtuSizing::clear();
    cellTypes_.clear();
    vertLabels_.clear();
    vertOffset_.clear();
    faceLabels_.clear();
    faceOffset_.clear();

    maps_.clear();
}


void Foam::vtk::vtuCells::repopulate(const polyMesh& mesh)
{
    // vtuSizing::reset() called prior to this method

    cellTypes_.resize(nFieldCells());
    vertLabels_.resize(sizeOf(output_, slotType::CELLS));
    vertOffset_.resize(sizeOf(output_, slotType::CELLS_OFFSETS));
    faceLabels_.resize(sizeOf(output_, slotType::FACES));
    faceOffset_.resize(sizeOf(output_, slotType::FACES_OFFSETS));

    switch (output_)
    {
        case contentType::LEGACY:
        {
            populateLegacy
            (
                mesh,
                cellTypes_,
                vertLabels_,
                maps_
            );
            break;
        }

        case contentType::XML:
        {
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
        }

        case contentType::INTERNAL1:
        case contentType::INTERNAL2:
        {
            populateInternal
            (
                mesh,
                cellTypes_,
                vertLabels_,
                vertOffset_,
                faceLabels_,
                faceOffset_,
                maps_,
                output_
            );
            break;
        }
    }
}


void Foam::vtk::vtuCells::reset(const polyMesh& mesh)
{
    vtuSizing::reset(mesh, decomposeRequest_);
    repopulate(mesh);
}


void Foam::vtk::vtuCells::reset
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


void Foam::vtk::vtuCells::renumberCells(const labelUList& mapping)
{
    maps_.renumberCells(mapping);
}


void Foam::vtk::vtuCells::renumberPoints(const labelUList& mapping)
{
    maps_.renumberPoints(mapping);
}


// ************************************************************************* //

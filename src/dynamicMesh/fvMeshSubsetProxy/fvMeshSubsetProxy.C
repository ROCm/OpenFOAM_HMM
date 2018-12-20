/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "fvMeshSubsetProxy.H"
#include "cellSet.H"
#include "cellZone.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshSubsetProxy::fvMeshSubsetProxy(fvMesh& baseMesh)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(-1),
    type_(NONE),
    name_(),
    names_()
{
    if (useSubMesh())
    {
        correct();
    }
}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    const subsetType type,
    const word& selectionName,
    label exposedPatchId
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(exposedPatchId),
    type_(selectionName.empty() ? NONE : type),
    name_(),
    names_()
{
    if (type_ == ZONES)
    {
        // Populate wordRes for ZONES
        names_.resize(1);
        names_.first() = selectionName;
    }
    else if (type_ != NONE)
    {
        name_ = selectionName;
    }

    if (useSubMesh())
    {
        correct();
    }
}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    const wordRes& zoneNames,
    label exposedPatchId
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(exposedPatchId),
    type_(ZONES),
    name_(),
    names_(zoneNames)
{
    if (useSubMesh())
    {
        correct();
    }
}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    wordRes&& zoneNames,
    label exposedPatchId
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(exposedPatchId),
    type_(ZONES),
    name_(),
    names_(std::move(zoneNames))
{
    if (useSubMesh())
    {
        correct();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshSubsetProxy::correct(bool verbose)
{
    if (type_ == NONE)
    {
        subsetter_.clear();
        return;
    }

    const label nCells = baseMesh_.nCells();

    bitSet selectedCells;

    if (type_ == SET)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellSet " << name_ << endl;
        }

        cellSet cset(baseMesh_, name_);

        selectedCells.resize(nCells);
        for (const label idx : cset)
        {
            selectedCells.set(idx);
        }
    }
    else if (type_ == ZONE)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellZone " << name_ << endl;
        }

        selectedCells.resize(nCells);
        selectedCells.set(baseMesh_.cellZones()[name_]);
    }
    else if (type_ == ZONES)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellZones "
                << flatOutput(names_) << endl;
        }

        selectedCells = baseMesh_.cellZones().selection(names_);
    }

    subsetter_.setCellSubset(selectedCells, exposedPatchId_);
}


Foam::polyMesh::readUpdateState Foam::fvMeshSubsetProxy::readUpdate()
{
    const polyMesh::readUpdateState meshState = baseMesh_.readUpdate();

    if
    (
        meshState == polyMesh::TOPO_CHANGE
     || meshState == polyMesh::TOPO_PATCH_CHANGE
    )
    {
        correct(true);
    }

    return meshState;
}


// ************************************************************************* //

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

#include "meshSubsetHelper.H"

#include "cellSet.H"
#include "cellZone.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshSubsetHelper::meshSubsetHelper
(
    fvMesh& baseMesh
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    type_(NONE),
    name_()
{
    correct();
}


Foam::meshSubsetHelper::meshSubsetHelper
(
    fvMesh& baseMesh,
    const subsetType type,
    const word& name
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    type_(name.empty() ? NONE : type),
    name_(name)
{
    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshSubsetHelper::correct(bool verbose)
{
    if (type_ == SET)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellSet " << name_ << endl;
        }

        cellSet subset(baseMesh_, name_);
        subsetter_.setLargeCellSubset(subset);
    }
    else if (type_ == ZONE)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellZone " << name_ << endl;
        }

        labelHashSet subset(baseMesh_.cellZones()[name_]);
        subsetter_.setLargeCellSubset(subset, 0);
    }
}


Foam::polyMesh::readUpdateState Foam::meshSubsetHelper::readUpdate()
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

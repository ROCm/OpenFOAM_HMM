/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "processorFaMeshes.H"
#include "Time.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::processorFaMeshes::read()
{
    forAll(fvMeshes_, proci)
    {
        meshes_.set(proci, new faMesh(fvMeshes_[proci]));

        // Read the addressing information

        IOobject ioAddr
        (
            "procAddressing",
            "constant",  // Placeholder
            faMesh::meshSubDir,
            meshes_[proci].thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        const auto& runTime = meshes_[proci].thisDb().time();
        const auto& meshDir = meshes_[proci].meshDir();

        // pointProcAddressing (faMesh)
        ioAddr.rename("pointProcAddressing");
        ioAddr.instance() = runTime.findInstance(meshDir, ioAddr.name());
        pointProcAddressing_.set(proci, new labelIOList(ioAddr));

        // edgeProcAddressing (faMesh)
        ioAddr.rename("edgeProcAddressing");
        ioAddr.instance() = runTime.findInstance(meshDir, ioAddr.name());
        edgeProcAddressing_.set(proci, new labelIOList(ioAddr));

        // faceProcAddressing (faMesh)
        ioAddr.rename("faceProcAddressing");
        ioAddr.instance() = runTime.findInstance(meshDir, ioAddr.name());
        faceProcAddressing_.set(proci, new labelIOList(ioAddr));

        // boundaryProcAddressing (faMesh)
        ioAddr.rename("boundaryProcAddressing");
        ioAddr.instance() = runTime.findInstance(meshDir, ioAddr.name());
        boundaryProcAddressing_.set(proci, new labelIOList(ioAddr));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorFaMeshes::processorFaMeshes
(
    const UPtrList<fvMesh>& procFvMeshes
)
:
    fvMeshes_(procFvMeshes),
    meshes_(procFvMeshes.size()),
    pointProcAddressing_(meshes_.size()),
    edgeProcAddressing_(meshes_.size()),
    faceProcAddressing_(meshes_.size()),
    boundaryProcAddressing_(meshes_.size())
{
    read();
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::processorFaMeshes::removeFiles(const faMesh& mesh)
{
    IOobject ioAddr
    (
        "procAddressing",
        mesh.facesInstance(),
        faMesh::meshSubDir,
        mesh.thisDb(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false  // not registered
    );

    // procAddressing
    rm(ioAddr.objectPath());

    // pointProcAddressing
    ioAddr.rename("pointProcAddressing");
    rm(ioAddr.objectPath());

    // edgeProcAddressing
    ioAddr.rename("edgeProcAddressing");
    rm(ioAddr.objectPath());

    // faceProcAddressing
    ioAddr.rename("faceProcAddressing");
    rm(ioAddr.objectPath());

    // boundaryProcAddressing
    ioAddr.rename("boundaryProcAddressing");
    rm(ioAddr.objectPath());
}


// ************************************************************************* //

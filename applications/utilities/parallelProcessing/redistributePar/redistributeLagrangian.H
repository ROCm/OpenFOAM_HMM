/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

Description
    Reading, reconstruct, redistribution of lagrangian fields.

\*---------------------------------------------------------------------------*/

#ifndef Foam_redistributeLagrangian_H
#define Foam_redistributeLagrangian_H

#include "parLagrangianDistributor.H"
#include "unmappedPassivePositionParticleCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Read clouds (note: might not be present on all processors)
PtrList<unmappedPassivePositionParticleCloud>
readLagrangian
(
    const fvMesh& mesh,
    const wordList& cloudNames,
    const wordRes& selectedFields
)
{
    PtrList<unmappedPassivePositionParticleCloud> clouds(cloudNames.size());

    if (!cloudNames.empty())
    {
        (void)mesh.tetBasePtIs();
    }

    // Setup clouds
    forAll(cloudNames, i)
    {
        //Pout<< "Loading cloud " << cloudNames[i] << endl;
        clouds.set
        (
            i,
            new unmappedPassivePositionParticleCloud(mesh, cloudNames[i], false)
        );

        //for (passivePositionParticle& p : clouds[i]))
        //{
        //    Pout<< "Particle position:" << p.position()
        //        << " cell:" << p.cell()
        //        << " with cc:" << mesh.cellCentres()[p.cell()]
        //        << endl;
        //}

        IOobjectList cloudObjs(clouds[i], clouds[i].time().timeName());

        //Pout<< "Found cloud objects:" << cloudObjs.names() << endl;

        parLagrangianDistributor::readAllFields
        (
            clouds[i],
            cloudObjs,
            selectedFields
        );
    }

    return clouds;
}


// Read clouds (note: might not be present on all processors)
PtrList<unmappedPassivePositionParticleCloud>
readLagrangian
(
    const fvMesh& mesh,
    const wordRes& selectedFields
)
{
    wordList cloudNames;
    List<wordList> fieldNames;
    // Find all cloudNames on all processors
    parLagrangianDistributor::findClouds(mesh, cloudNames, fieldNames);

    return readLagrangian(mesh, cloudNames, selectedFields);
}


void reconstructLagrangian
(
    autoPtr<parLagrangianDistributor>& distributorPtr,
    const fvMesh& baseMesh,
    const fvMesh& mesh,
    const mapDistributePolyMesh& distMap,
    const wordRes& selectedFields
)
{
    // Clouds (note: might not be present on all processors)

    wordList cloudNames;
    List<wordList> fieldNames;
    // Find all cloudNames on all processors
    parLagrangianDistributor::findClouds(mesh, cloudNames, fieldNames);

    if (cloudNames.empty())
    {
        // Nothing to do
        return;
    }

    // Use existing or create distributor
    if (!distributorPtr)
    {
        distributorPtr.reset
        (
            new parLagrangianDistributor
            (
                mesh,
                baseMesh,
                mesh.nCells(),      // range of cell indices in clouds
                distMap
            )
        );
    }
    const auto& distributor = *distributorPtr;

    for (const word& cloudName : cloudNames)
    {
        Info<< "Reconstructing lagrangian fields for cloud "
            << cloudName << nl << endl;

        IOobjectList cloudObjs
        (
            mesh,
            mesh.time().timeName(),
            cloud::prefix/cloudName
        );

        autoPtr<mapDistributeBase> lagrangianMapPtr =
            distributor.distributeLagrangianPositions
            (
                cloudName
            );

        distributor.distributeAllFields
        (
            lagrangianMapPtr(),
            cloudName,
            cloudObjs,
            selectedFields
        );
    }
}


void redistributeLagrangian
(
    autoPtr<parLagrangianDistributor>& distributorPtr,
    const fvMesh& mesh,
    const label nOldCells,
    const mapDistributePolyMesh& distMap,
    PtrList<unmappedPassivePositionParticleCloud>& clouds
)
{
    if (clouds.empty())
    {
        // Nothing to do
        return;
    }

    // Use existing or create distributor
    if (!distributorPtr)
    {
        distributorPtr.reset
        (
            new parLagrangianDistributor
            (
                mesh,
                mesh,
                nOldCells,  // range of cell indices in clouds
                distMap
            )
        );
    }
    const auto& distributor = *distributorPtr;

    for (auto& cloud : clouds)
    {
        autoPtr<mapDistributeBase> lagrangianMapPtr =
            distributor.distributeLagrangianPositions(cloud);

        distributor.distributeAllStoredFields
        (
            lagrangianMapPtr(),
            cloud
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //

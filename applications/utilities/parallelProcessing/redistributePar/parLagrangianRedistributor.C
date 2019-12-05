/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "ListOps.H"
#include "parLagrangianRedistributor.H"
#include "passivePositionParticleCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parLagrangianRedistributor::parLagrangianRedistributor
(
    const fvMesh& srcMesh,
    const fvMesh& tgtMesh,
    const label nSrcCells,
    const mapDistributePolyMesh& distMap
)
:
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh),
    distMap_(distMap)
{
    const mapDistribute& cellMap = distMap_.cellMap();

    // Get destination processors and cells
    destinationProcID_ = labelList(tgtMesh_.nCells(), Pstream::myProcNo());
    cellMap.reverseDistribute(nSrcCells, destinationProcID_);

    destinationCell_ = identity(tgtMesh_.nCells());
    cellMap.reverseDistribute(nSrcCells, destinationCell_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Find all clouds (on all processors) and for each cloud all the objects.
// Result will be synchronised on all processors
void Foam::parLagrangianRedistributor::findClouds
(
    const fvMesh& mesh,
    wordList& cloudNames,
    List<wordList>& objectNames
)
{
    fileNameList localCloudDirs
    (
        readDir
        (
            mesh.time().timePath()
          / mesh.dbDir()
          / cloud::prefix,
            fileName::DIRECTORY
        )
    );

    cloudNames.setSize(localCloudDirs.size());
    forAll(localCloudDirs, i)
    {
        cloudNames[i] = localCloudDirs[i];
    }

    // Synchronise cloud names
    Pstream::combineGather(cloudNames, ListOps::uniqueEqOp<word>());
    Pstream::combineScatter(cloudNames);

    objectNames.setSize(cloudNames.size());

    for (const fileName& localCloudName : localCloudDirs)
    {
        // Do local scan for valid cloud objects
        IOobjectList localObjs
        (
            mesh,
            mesh.time().timeName(),
            cloud::prefix/localCloudName
        );

        bool isCloud = false;
        if (localObjs.erase("coordinates"))
        {
            isCloud = true;
        }
        if (localObjs.erase("positions"))
        {
            isCloud = true;
        }

        if (isCloud)
        {
            // Has coordinates/positions - so must be a valid cloud

            const label cloudi = cloudNames.find(localCloudName);

            objectNames[cloudi] = localObjs.sortedNames();
        }
    }

    // Synchronise objectNames
    forAll(objectNames, i)
    {
        Pstream::combineGather(objectNames[i], ListOps::uniqueEqOp<word>());
        Pstream::combineScatter(objectNames[i]);
    }
}


Foam::autoPtr<Foam::mapDistributeBase>
Foam::parLagrangianRedistributor::redistributeLagrangianPositions
(
    passivePositionParticleCloud& lpi
) const
{
    //Debug(lpi.size());

    const label oldLpi = lpi.size();

    labelListList subMap;


    // Allocate transfer buffers
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    {
        // List of lists of particles to be transferred for all of the
        // neighbour processors
        List<IDLList<passivePositionParticle>> particleTransferLists
        (
            Pstream::nProcs()
        );

        // Per particle the destination processor
        labelList destProc(lpi.size());

        label particleI = 0;
        for (passivePositionParticle& ppi : lpi)
        {
            const label destProcI = destinationProcID_[ppi.cell()];
            const label destCellI = destinationCell_[ppi.cell()];

            ppi.cell() = destCellI;
            destProc[particleI++] = destProcI;
            particleTransferLists[destProcI].append(lpi.remove(&ppi));
        }


        // Per processor the indices of the particles to send
        subMap = invertOneToMany(Pstream::nProcs(), destProc);


        // Stream into send buffers
        forAll(particleTransferLists, procI)
        {
            //Pout<< "To proc " << procI << " sending "
            //    << particleTransferLists[procI] << endl;
            if (particleTransferLists[procI].size())
            {
                UOPstream particleStream(procI, pBufs);
                particleStream << particleTransferLists[procI];
            }
        }
    }


    // Start sending. Sets number of bytes transferred
    labelList allNTrans(Pstream::nProcs());
    pBufs.finishedSends(allNTrans);


    {
        // Temporarily rename original cloud so we can construct a new one
        // (to distribute the positions) without getting a duplicate
        // registration warning
        const word cloudName = lpi.name();
        lpi.rename(cloudName + "_old");

        // New cloud on tgtMesh
        passivePositionParticleCloud lagrangianPositions
        (
            tgtMesh_,
            cloudName,
            IDLList<passivePositionParticle>()
        );


        // Retrieve from receive buffers
        forAll(allNTrans, procI)
        {
            const label nRec = allNTrans[procI];

            //Pout<< "From processor " << procI << " receiving bytes " << nRec
            //    << endl;

            if (nRec)
            {
                UIPstream particleStream(procI, pBufs);

                // Receive particles and locate them
                IDLList<passivePositionParticle> newParticles
                (
                    particleStream,
                    passivePositionParticle::iNew(tgtMesh_)
                );

                for (passivePositionParticle& newp : newParticles)
                {
                    lagrangianPositions.addParticle(newParticles.remove(&newp));
                }
            }
        }

        if (lagrangianPositions.size())
        {
            // Write coordinates file
            IOPosition<passivePositionParticleCloud>
            (
                lagrangianPositions
            ).write();

            // Optionally write positions file in v1706 format and earlier
            if (particle::writeLagrangianPositions)
            {
                IOPosition<passivePositionParticleCloud>
                (
                     lagrangianPositions,
                     cloud::geometryType::POSITIONS
                 ).write();
            }
        }
        else if (oldLpi)
        {
            // When running with -overwrite it should also delete the old
            // files. Below works but is not optimal.

            // Remove any existing coordinates
            const fileName oldCoords
            (
                IOPosition<passivePositionParticleCloud>
                (
                    lagrangianPositions
                ).objectPath()
            );
            Foam::rm(oldCoords);

            // Remove any existing positions
            const fileName oldPos
            (
                IOPosition<passivePositionParticleCloud>
                (
                    lagrangianPositions,
                    cloud::geometryType::POSITIONS
                ).objectPath()
            );
            Foam::rm(oldPos);
        }

        // Restore cloud name
        lpi.rename(cloudName);
    }

    // Work the send indices (subMap) into a mapDistributeBase
    labelListList sizes(Pstream::nProcs());
    labelList& nsTransPs = sizes[Pstream::myProcNo()];
    nsTransPs.setSize(Pstream::nProcs());
    forAll(subMap, sendProcI)
    {
        nsTransPs[sendProcI] = subMap[sendProcI].size();
    }
    // Send sizes across. Note: blocks.
    combineReduce(sizes, Pstream::listEq());

    labelListList constructMap(Pstream::nProcs());
    label constructSize = 0;
    forAll(constructMap, procI)
    {
        const label nRecv = sizes[procI][UPstream::myProcNo()];

        labelList& map = constructMap[procI];

        map.setSize(nRecv);
        forAll(map, i)
        {
            map[i] = constructSize++;
        }
    }


    return autoPtr<mapDistributeBase>::New
    (
        constructSize,
        std::move(subMap),
        std::move(constructMap)
    );
}


Foam::autoPtr<Foam::mapDistributeBase>
Foam::parLagrangianRedistributor::redistributeLagrangianPositions
(
    const word& cloudName
) const
{
    // Load cloud and send particle
    passivePositionParticleCloud lpi(srcMesh_, cloudName, false);

    return redistributeLagrangianPositions(lpi);
}


// ************************************************************************* //

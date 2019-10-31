/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

Author
    Hrvoje Jasak

\*----------------------------------------------------------------------------*/

#include "faGlobalMeshData.H"
#include "faMesh.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faGlobalMeshData::faGlobalMeshData(const faMesh& mesh)
:
    faProcessorTopology(mesh.boundary(), UPstream::worldComm),
    mesh_(mesh),
    nGlobalPoints_(-1),
    sharedPointLabels_(0),
    sharedPointAddr_(0)
{
    updateMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faGlobalMeshData::~faGlobalMeshData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::faMesh& Foam::faGlobalMeshData::mesh() const
{
    return mesh_;
}


void Foam::faGlobalMeshData::updateMesh()
{
    label polyMeshNGlobalPoints = mesh_().globalData().nGlobalPoints();

    const labelList& polyMeshSharedPointLabels =
        mesh_().globalData().sharedPointLabels();

    const labelList& polyMeshSharedPointAddr =
        mesh_().globalData().sharedPointAddr();

    labelHashSet sharedPointLabels;

    labelField globalList(polyMeshNGlobalPoints, Zero);

    forAll(mesh_.boundary(), patchI)
    {
        const faPatch& fap = mesh_.boundary()[patchI];

        if (isA<processorFaPatch>(fap))
        {
            const labelList& localPointLabels = fap.pointLabels();

            forAll(localPointLabels, pointI)
            {
                label polyMeshPoint =
                    mesh_.patch().meshPoints()[localPointLabels[pointI]];

                const label sharedPolyMeshPoint =
                    polyMeshSharedPointLabels.find(polyMeshPoint);

                if
                (
                    sharedPolyMeshPoint != -1
                 && !sharedPointLabels.found(localPointLabels[pointI])
                )
                {
                    globalList[polyMeshSharedPointAddr[sharedPolyMeshPoint]]
                        += 1;

                    sharedPointLabels.insert(localPointLabels[pointI]);
                }
            }
        }
    }

    sharedPointLabels_ = sharedPointLabels.toc();

    combineReduce(globalList, plusEqOp<labelField>());

    nGlobalPoints_ = 0;
    for (label i=0; i<globalList.size(); ++i)
    {
        if (globalList[i] > 0)
        {
            globalList[i] = ++nGlobalPoints_;
        }
    }

    sharedPointAddr_.setSize(sharedPointLabels_.size());
    forAll(sharedPointAddr_, pointI)
    {
        const label polyMeshSharedPointIndex =
            polyMeshSharedPointLabels.find
            (
                mesh_.patch().meshPoints()[sharedPointLabels_[pointI]]
            );

        sharedPointAddr_[pointI] =
            globalList[polyMeshSharedPointAddr[polyMeshSharedPointIndex]]
          - 1;
    }
}


// ************************************************************************* //

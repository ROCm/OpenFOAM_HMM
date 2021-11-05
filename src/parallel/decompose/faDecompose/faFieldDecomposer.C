/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "faFieldDecomposer.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const label sizeBeforeMapping,
    const labelUList& addressingSlice,
    const label addressingOffset
)
:
    sizeBeforeMapping_(sizeBeforeMapping),
    directAddressing_(addressingSlice)
{
    forAll(directAddressing_, i)
    {
        // Subtract one to align addressing.
        // directAddressing_[i] -= addressingOffset + 1;
        // ZT, 12/Nov/2010
        directAddressing_[i] -= addressingOffset;
    }
}


Foam::faFieldDecomposer::processorAreaPatchFieldDecomposer::
processorAreaPatchFieldDecomposer
(
    const label nTotalFaces,
    const labelUList& owner,  // == mesh.edgeOwner()
    const labelUList& neigh,  // == mesh.edgeNeighbour()
    const labelUList& addressingSlice,
    const scalarField& weights
)
:
    sizeBeforeMapping_(nTotalFaces),
    addressing_(addressingSlice.size()),
    weights_(addressingSlice.size())
{
    forAll(addressing_, i)
    {
        // Subtract one to align addressing.
        label ai = addressingSlice[i];
//         label ai = mag(addressingSlice[i]) - 1;

        if (ai < neigh.size())
        {
            // This is a regular edge. it has been an internal edge
            // of the original mesh and now it has become a edge
            // on the parallel boundary
            addressing_[i].resize(2);
            weights_[i].resize(2);

            addressing_[i][0] = owner[ai];
            addressing_[i][1] = neigh[ai];

            if (ai < weights.size())
            {
                // Edge weights exist/are usable
                weights_[i][0] = weights[ai];
                weights_[i][1] = 1.0 - weights[ai];
            }
            else
            {
                // No edge weights. use equal weighting
                weights_[i][0] = 0.5;
                weights_[i][1] = 0.5;
            }
        }
        else
        {
            // This is a edge that used to be on a cyclic boundary
            // but has now become a parallel patch edge. I cannot
            // do the interpolation properly (I would need to look
            // up the different (edge) list of data), so I will
            // just grab the value from the owner face

            addressing_[i].resize(1);
            weights_[i].resize(1);

            addressing_[i][0] = owner[ai];

            weights_[i][0] = 1.0;
        }
    }
}


Foam::faFieldDecomposer::processorAreaPatchFieldDecomposer::
processorAreaPatchFieldDecomposer
(
    const faMesh& mesh,
    const labelUList& addressingSlice
)
:
    processorAreaPatchFieldDecomposer
    (
        mesh.nFaces(),
        mesh.edgeOwner(),
        mesh.edgeNeighbour(),
        addressingSlice,
        mesh.weights().internalField()
    )
{}


Foam::faFieldDecomposer::processorEdgePatchFieldDecomposer::
processorEdgePatchFieldDecomposer
(
    label sizeBeforeMapping,
    const labelUList& addressingSlice
)
:
    sizeBeforeMapping_(sizeBeforeMapping),
    addressing_(addressingSlice.size()),
    weights_(addressingSlice.size())
{
    forAll(addressing_, i)
    {
        addressing_[i].resize(1);
        weights_[i].resize(1);

        addressing_[i][0] = mag(addressingSlice[i]) - 1;
        weights_[i][0] = sign(addressingSlice[i]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faFieldDecomposer::faFieldDecomposer
(
    const Foam::zero,
    const faMesh& procMesh,
    const labelList& edgeAddressing,
    const labelList& faceAddressing,
    const labelList& boundaryAddressing
)
:
    procMesh_(procMesh),
    edgeAddressing_(edgeAddressing),
    faceAddressing_(faceAddressing),
    boundaryAddressing_(boundaryAddressing),
    // Mappers
    patchFieldDecomposerPtrs_(),
    processorAreaPatchFieldDecomposerPtrs_(),
    processorEdgePatchFieldDecomposerPtrs_()
{}


Foam::faFieldDecomposer::faFieldDecomposer
(
    const faMesh& completeMesh,
    const faMesh& procMesh,
    const labelList& edgeAddressing,
    const labelList& faceAddressing,
    const labelList& boundaryAddressing
)
:
    faFieldDecomposer
    (
        zero{},
        procMesh,
        edgeAddressing,
        faceAddressing,
        boundaryAddressing
    )
{
    reset(completeMesh);
}


Foam::faFieldDecomposer::faFieldDecomposer
(
    const label nTotalFaces,
    const List<labelRange>& boundaryRanges,
    const labelUList& edgeOwner,
    const labelUList& edgeNeigbour,

    const faMesh& procMesh,
    const labelList& edgeAddressing,
    const labelList& faceAddressing,
    const labelList& boundaryAddressing
)
:
    faFieldDecomposer
    (
        zero{},
        procMesh,
        edgeAddressing,
        faceAddressing,
        boundaryAddressing
    )
{
    reset(nTotalFaces, boundaryRanges, edgeOwner, edgeNeigbour);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faFieldDecomposer::empty() const
{
    return patchFieldDecomposerPtrs_.empty();
}


void Foam::faFieldDecomposer::clear()
{
    patchFieldDecomposerPtrs_.clear();
    processorAreaPatchFieldDecomposerPtrs_.clear();
    processorEdgePatchFieldDecomposerPtrs_.clear();
}


void Foam::faFieldDecomposer::reset
(
    const label nTotalFaces,
    const List<labelRange>& boundaryRanges,
    const labelUList& edgeOwner,
    const labelUList& edgeNeigbour
)
{
    clear();
    const label nMappers = procMesh_.boundary().size();

    patchFieldDecomposerPtrs_.resize(nMappers);
    processorAreaPatchFieldDecomposerPtrs_.resize(nMappers);
    processorEdgePatchFieldDecomposerPtrs_.resize(nMappers);

    forAll(boundaryAddressing_, patchi)
    {
        const label oldPatchi = boundaryAddressing_[patchi];
        const faPatch& fap = procMesh_.boundary()[patchi];
        const labelSubList localPatchSlice(fap.patchSlice(edgeAddressing_));

        if (oldPatchi >= 0)
        {
            patchFieldDecomposerPtrs_.set
            (
                patchi,
                new patchFieldDecomposer
                (
                    boundaryRanges[oldPatchi].size(),
                    localPatchSlice,
                    boundaryRanges[oldPatchi].start()
                )
            );
        }
        else
        {
            processorAreaPatchFieldDecomposerPtrs_.set
            (
                patchi,
                new processorAreaPatchFieldDecomposer
                (
                    nTotalFaces,
                    edgeOwner,
                    edgeNeigbour,
                    localPatchSlice
                )
            );

            processorEdgePatchFieldDecomposerPtrs_.set
            (
                patchi,
                new processorEdgePatchFieldDecomposer
                (
                    procMesh_.boundary()[patchi].size(),
                    localPatchSlice
                )
            );
        }
    }
}


void Foam::faFieldDecomposer::reset(const faMesh& completeMesh)
{
    clear();
    const label nMappers = procMesh_.boundary().size();
    patchFieldDecomposerPtrs_.resize(nMappers);
    processorAreaPatchFieldDecomposerPtrs_.resize(nMappers);
    processorEdgePatchFieldDecomposerPtrs_.resize(nMappers);

    // Create weightings now - needed for proper parallel synchronization
    (void)completeMesh.weights();

    // faPatches don't have their own start() - so these are invariant
    const labelList completePatchStarts
    (
        completeMesh.boundary().patchStarts()
    );

    forAll(boundaryAddressing_, patchi)
    {
        const label oldPatchi = boundaryAddressing_[patchi];
        const faPatch& fap = procMesh_.boundary()[patchi];
        const labelSubList localPatchSlice(fap.patchSlice(edgeAddressing_));

        if (oldPatchi >= 0)
        {
            patchFieldDecomposerPtrs_.set
            (
                patchi,
                new patchFieldDecomposer
                (
                    completeMesh.boundary()[oldPatchi].size(),
                    localPatchSlice,
                    completePatchStarts[oldPatchi]
                )
            );
        }
        else
        {
            processorAreaPatchFieldDecomposerPtrs_.set
            (
                patchi,
                new processorAreaPatchFieldDecomposer
                (
                    completeMesh,
                    localPatchSlice
                )
            );

            processorEdgePatchFieldDecomposerPtrs_.set
            (
                patchi,
                new processorEdgePatchFieldDecomposer
                (
                    procMesh_.boundary()[patchi].size(),
                    localPatchSlice
                )
            );
        }
    }
}


// ************************************************************************* //

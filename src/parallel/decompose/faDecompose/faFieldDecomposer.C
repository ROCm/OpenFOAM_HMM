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
    const faMesh& mesh,
    const labelUList& addressingSlice
)
:
    sizeBeforeMapping_(mesh.nFaces()),
    addressing_(addressingSlice.size()),
    weights_(addressingSlice.size())
{
    const scalarField& weights = mesh.weights().internalField();
    const labelList& own = mesh.edgeOwner();
    const labelList& neighb = mesh.edgeNeighbour();

    forAll(addressing_, i)
    {
        // Subtract one to align addressing.
        label ai = addressingSlice[i];
//         label ai = mag(addressingSlice[i]) - 1;

        if (ai < neighb.size())
        {
            // This is a regular edge. it has been an internal edge
            // of the original mesh and now it has become a edge
            // on the parallel boundary
            addressing_[i].setSize(2);
            weights_[i].setSize(2);

            addressing_[i][0] = own[ai];
            addressing_[i][1] = neighb[ai];

            weights_[i][0] = weights[ai];
            weights_[i][1] = 1.0 - weights[ai];
        }
        else
        {
            // This is a edge that used to be on a cyclic boundary
            // but has now become a parallel patch edge. I cannot
            // do the interpolation properly (I would need to look
            // up the different (edge) list of data), so I will
            // just grab the value from the owner face
            //
            addressing_[i].setSize(1);
            weights_[i].setSize(1);

            addressing_[i][0] = own[ai];

            weights_[i][0] = 1.0;
        }
    }
}


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
        addressing_[i].setSize(1);
        weights_[i].setSize(1);

        addressing_[i][0] = mag(addressingSlice[i]) - 1;
        weights_[i][0] = sign(addressingSlice[i]);
    }
}


Foam::faFieldDecomposer::faFieldDecomposer
(
    const faMesh& completeMesh,
    const faMesh& procMesh,
    const labelList& edgeAddressing,
    const labelList& faceAddressing,
    const labelList& boundaryAddressing
)
:
    completeMesh_(completeMesh),
    procMesh_(procMesh),
    edgeAddressing_(edgeAddressing),
    faceAddressing_(faceAddressing),
    boundaryAddressing_(boundaryAddressing),

    patchFieldDecomposerPtrs_(procMesh_.boundary().size()),
    processorAreaPatchFieldDecomposerPtrs_(procMesh_.boundary().size()),
    processorEdgePatchFieldDecomposerPtrs_(procMesh_.boundary().size())
{
    forAll(boundaryAddressing_, patchi)
    {
        const label oldPatchi = boundaryAddressing_[patchi];

        if (oldPatchi >= 0)
        {
            patchFieldDecomposerPtrs_.set
            (
                patchi,
                new patchFieldDecomposer
                (
                    completeMesh_.boundary()[oldPatchi].size(),
                    procMesh_.boundary()[patchi].patchSlice(edgeAddressing_),
                    completeMesh_.boundary()[oldPatchi].start()
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
                    completeMesh_,
                    procMesh_.boundary()[patchi].patchSlice(edgeAddressing_)
                )
            );

            processorEdgePatchFieldDecomposerPtrs_.set
            (
                patchi,
                new processorEdgePatchFieldDecomposer
                (
                    procMesh_.boundary()[patchi].size(),
                    static_cast<const labelUList&>
                    (
                        procMesh_.boundary()[patchi].patchSlice
                        (
                            edgeAddressing_
                        )
                    )
                )
            );
        }
    }
}


// ************************************************************************* //

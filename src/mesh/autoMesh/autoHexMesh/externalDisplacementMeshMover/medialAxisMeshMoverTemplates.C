/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "medialAxisMeshMover.H"
#include "syncTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::medialAxisMeshMover::weightedSum
(
    const polyMesh& mesh,
    const PackedBoolList& isMasterEdge,
    const labelList& meshEdges,
    const labelList& meshPoints,
    const edgeList& edges,
    const scalarField& edgeWeights,
    const Field<Type>& pointData,
    Field<Type>& sum
)
{
    if
    (
        mesh.nEdges() != isMasterEdge.size()
     || edges.size() != meshEdges.size()
     || edges.size() != edgeWeights.size()
     || meshPoints.size() != pointData.size()
    )
    {
        FatalErrorIn("medialAxisMeshMover::weightedSum(..)")
            << "Inconsistent sizes for edge or point data:"
            << " meshEdges:" << meshEdges.size()
            << " isMasterEdge:" << isMasterEdge.size()
            << " edgeWeights:" << edgeWeights.size()
            << " edges:" << edges.size()
            << " pointData:" << pointData.size()
            << " meshPoints:" << meshPoints.size()
            << abort(FatalError);
    }

    sum.setSize(meshPoints.size());
    sum = pTraits<Type>::zero;

    forAll(edges, edgeI)
    {
        if (isMasterEdge.get(meshEdges[edgeI]) == 1)
        {
            const edge& e = edges[edgeI];

            scalar eWeight = edgeWeights[edgeI];

            label v0 = e[0];
            label v1 = e[1];

            sum[v0] += eWeight*pointData[v1];
            sum[v1] += eWeight*pointData[v0];
        }
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        sum,
        plusEqOp<Type>(),
        pTraits<Type>::zero     // null value
    );
}


// ************************************************************************* //

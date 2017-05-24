/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "edgeMesh.H"
#include "mergePoints.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"
#include "ListOps.H"
#include "EdgeMap.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(edgeMesh, 0);
    defineRunTimeSelectionTable(edgeMesh, fileExtension);
    defineMemberFunctionSelectionTable(edgeMesh,write,fileExtension);
}


Foam::wordHashSet Foam::edgeMesh::readTypes()
{
    return wordHashSet(*fileExtensionConstructorTablePtr_);
}


Foam::wordHashSet Foam::edgeMesh::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::edgeMesh::canReadType
(
    const word& ext,
    const bool verbose
)
{
    return checkSupport
    (
        readTypes(),
        ext,
        verbose,
        "reading"
   );
}


bool Foam::edgeMesh::canWriteType
(
    const word& ext,
    const bool verbose
)
{
    return checkSupport
    (
        writeTypes(),
        ext,
        verbose,
        "writing"
    );
}


bool Foam::edgeMesh::canRead
(
    const fileName& name,
    const bool verbose
)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }
    return canReadType(ext, verbose);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::edgeMesh::calcPointEdges() const
{
    if (pointEdgesPtr_.valid())
    {
        FatalErrorInFunction
            << "pointEdges already calculated." << abort(FatalError);
    }

    pointEdgesPtr_.reset(new labelListList(points_.size()));
    labelListList& pointEdges = pointEdgesPtr_();

    invertManyToMany(pointEdges.size(), edges_, pointEdges);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeMesh::edgeMesh()
:
    fileFormats::edgeMeshFormatsCore(),
    points_(0),
    edges_(0),
    pointEdgesPtr_(nullptr)
{}


Foam::edgeMesh::edgeMesh
(
    const pointField& points,
    const edgeList& edges
)
:
    fileFormats::edgeMeshFormatsCore(),
    points_(points),
    edges_(edges),
    pointEdgesPtr_(nullptr)
{}


Foam::edgeMesh::edgeMesh
(
    const Xfer<pointField>& pointLst,
    const Xfer<edgeList>& edgeLst
)
:
    fileFormats::edgeMeshFormatsCore(),
    points_(0),
    edges_(0),
    pointEdgesPtr_(nullptr)
{
    points_.transfer(pointLst());
    edges_.transfer(edgeLst());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::edgeMesh::~edgeMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::edgeMesh::clear()
{
    points_.clear();
    edges_.clear();
    pointEdgesPtr_.clear();
}


void Foam::edgeMesh::reset
(
    const Xfer<pointField>& pointLst,
    const Xfer<edgeList>& edgeLst
)
{
    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (notNull(pointLst))
    {
        points_.transfer(pointLst());
    }

    if (notNull(edgeLst))
    {
        edges_.transfer(edgeLst());

        // connectivity likely changed
        pointEdgesPtr_.clear();
    }
}


void Foam::edgeMesh::transfer(edgeMesh& mesh)
{
    points_.transfer(mesh.points_);
    edges_.transfer(mesh.edges_);
    pointEdgesPtr_ = mesh.pointEdgesPtr_;
}


Foam::Xfer<Foam::edgeMesh> Foam::edgeMesh::xfer()
{
    return xferMove(*this);
}


Foam::label Foam::edgeMesh::regions(labelList& edgeRegion) const
{
    edgeRegion.setSize(edges_.size());
    edgeRegion = -1;

    label startEdgeI = 0;
    label currentRegion = 0;

    while (true)
    {
        while (startEdgeI < edges_.size() && edgeRegion[startEdgeI] != -1)
        {
            startEdgeI++;
        }

        if (startEdgeI == edges_.size())
        {
            break;
        }

        // Found edge that has not yet been assigned a region.
        // Mark connected region with currentRegion starting at startEdgeI.

        edgeRegion[startEdgeI] = currentRegion;
        labelList edgesToVisit(1, startEdgeI);

        while (edgesToVisit.size())
        {
            // neighbours of current edgesToVisit
            DynamicList<label> newEdgesToVisit(edgesToVisit.size());

            // Mark all point connected edges with current region.
            forAll(edgesToVisit, i)
            {
                label edgeI = edgesToVisit[i];

                // Mark connected edges
                const edge& e = edges_[edgeI];

                forAll(e, fp)
                {
                    const labelList& pEdges = pointEdges()[e[fp]];

                    forAll(pEdges, pEdgeI)
                    {
                        label nbrEdgeI = pEdges[pEdgeI];

                        if (edgeRegion[nbrEdgeI] == -1)
                        {
                            edgeRegion[nbrEdgeI] = currentRegion;
                            newEdgesToVisit.append(nbrEdgeI);
                        }
                    }
                }
            }

            edgesToVisit.transfer(newEdgesToVisit);
        }

        currentRegion++;
    }
    return currentRegion;
}


void Foam::edgeMesh::scalePoints(const scalar scaleFactor)
{
    // avoid bad scaling
    if (scaleFactor > 0 && scaleFactor != 1.0)
    {
        points_ *= scaleFactor;
    }
}


void Foam::edgeMesh::mergePoints(const scalar mergeDist)
{
    pointField newPoints;
    labelList pointMap;

    const bool hasMerged = Foam::mergePoints
    (
        points_,
        mergeDist,
        false,
        pointMap,
        newPoints,
        vector::zero
    );

    if (hasMerged)
    {
        pointEdgesPtr_.clear();   // connectivity change

        points_.transfer(newPoints);

        forAll(edges_, edgeI)
        {
            edge& e = edges_[edgeI];

            e[0] = pointMap[e[0]];
            e[1] = pointMap[e[1]];
        }
    }

    this->mergeEdges();
}


void Foam::edgeMesh::mergeEdges()
{
    HashSet<edge, Hash<edge>> uniqEdges(2*edges_.size());
    PackedBoolList pointIsUsed(points_.size());

    label nUniqEdges = 0;
    label nUniqPoints = 0;
    forAll(edges_, edgeI)
    {
        const edge& e = edges_[edgeI];

        // Remove degenerate and repeated edges
        // - reordering (e[0] < e[1]) is not really necessary
        if (e[0] != e[1] && uniqEdges.insert(e))
        {
            if (nUniqEdges != edgeI)
            {
                edges_[nUniqEdges] = e;
            }
            edges_[nUniqEdges].sort();
            ++nUniqEdges;

            if (pointIsUsed.set(e[0], 1))
            {
                ++nUniqPoints;
            }
            if (pointIsUsed.set(e[1], 1))
            {
                ++nUniqPoints;
            }
        }
    }

    if (debug)
    {
        Info<< "Merging duplicate edges: "
            << (edges_.size() - nUniqEdges)
            << " edges will be deleted, "
            << (points_.size() - nUniqPoints)
            << " unused points will be removed." << endl;
    }

    if (nUniqEdges < edges_.size())
    {
        pointEdgesPtr_.clear(); // connectivity change
        edges_.setSize(nUniqEdges);  // truncate
    }

    if (nUniqPoints < points_.size())
    {
        pointEdgesPtr_.clear(); // connectivity change

        // build a oldToNew point-map and rewrite the points.
        // We can do this simultaneously since the point order is unchanged
        // and we are only effectively eliminating some entries.
        labelList pointMap(points_.size(), -1);

        label newId = 0;
        forAll(pointMap, pointi)
        {
            if (pointIsUsed[pointi])
            {
                pointMap[pointi] = newId;

                if (newId < pointi)
                {
                    // copy down
                    points_[newId] = points_[pointi];
                }
                ++newId;
            }
        }
        points_.setSize(newId);

        // Renumber edges - already sorted (above)
        forAll(edges_, edgeI)
        {
            edge& e = edges_[edgeI];

            e[0] = pointMap[e[0]];
            e[1] = pointMap[e[1]];
        }
    }

}


// ************************************************************************* //

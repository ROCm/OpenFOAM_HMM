/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "triSurfaceMesh.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(triSurfaceMesh, 0);
addToRunTimeSelectionTable(searchableSurface, triSurfaceMesh, dict);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Check file existence
const Foam::fileName& Foam::triSurfaceMesh::checkFile
(
    const fileName& fName,
    const fileName& objectName
)
{
    if (fName == fileName::null)
    {
        FatalErrorIn
        (
            "triSurfaceMesh::checkFile(const fileName&, const fileName&)"
        )   << "Cannot find triSurfaceMesh starting from "
            << objectName << exit(FatalError);
    }
    return fName;
}


const Foam::indexedOctree<Foam::treeDataTriSurface>&
    Foam::triSurfaceMesh::tree() const
{
    if (!tree_.valid())
    {
        treeBoundBox bb(points(), meshPoints());

        // Random number generator. Bit dodgy since not exactly random ;-)
        Random rndGen(65431);

        tree_.reset
        (
            new indexedOctree<treeDataTriSurface>
            (
                treeDataTriSurface(*this),
                bb.extend(rndGen, 1E-3),    // slightly randomize bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    return tree_();
}


const Foam::indexedOctree<Foam::treeDataEdge>&
    Foam::triSurfaceMesh::edgeTree() const
{
    if (!edgeTree_.valid())
    {
        treeBoundBox bb(localPoints());

        // Boundary edges
        labelList bEdges
        (
            identity
            (
                nEdges()
               -nInternalEdges()
            )
          + nInternalEdges()
        );

        // Random number generator. Bit dodgy since not exactly random ;-)
        Random rndGen(65431);

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    false,          // cachebb
                    edges(),        // edges
                    localPoints(),  // points
                    bEdges          // selected edges
                ),
                bb.extend(rndGen, 1E-3),    // slightly randomize bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }
    return edgeTree_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from triangles, patches, points.
Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const triSurface& s)
:
    searchableSurface(io.name()),
    objectRegistry(io),
    triSurface(s)
{}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io)
:
    searchableSurface(io.name()),
    objectRegistry(io),
    triSurface(checkFile(filePath(), objectPath()))
{}


Foam::triSurfaceMesh::triSurfaceMesh
(
    const word& name,
    const objectRegistry& obj,
    const dictionary& dict
)
:
    searchableSurface(name),
    objectRegistry
    (
        IOobject
        (
            name,
            "constant",
            "triSurface",
            obj,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    triSurface(checkFile(filePath(), objectPath()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::~triSurfaceMesh()
{}


void Foam::triSurfaceMesh::clearOut()
{
    tree_.clear();
    edgeTree_.clear();
    triSurface::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurfaceMesh::movePoints(const pointField& newPoints)
{
    tree_.clear();
    edgeTree_.clear();
    triSurface::movePoints(newPoints);    
}


Foam::pointIndexHit Foam::triSurfaceMesh::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    return tree().findNearest(sample, nearestDistSqr);
}


Foam::pointIndexHit Foam::triSurfaceMesh::findNearestOnEdge
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    return edgeTree().findNearest(sample, nearestDistSqr);
}


Foam::pointIndexHit Foam::triSurfaceMesh::findNearest
(
    const linePointRef& ln,
    treeBoundBox& tightest,
    point& linePoint
) const
{
    return tree().findNearest(ln, tightest, linePoint);
}


Foam::pointIndexHit Foam::triSurfaceMesh::findLine
(
    const point& start,
    const point& end
) const
{
    return tree().findLine(start, end);
}


Foam::pointIndexHit Foam::triSurfaceMesh::findLineAny
(
    const point& start,
    const point& end
) const
{
    return tree().findLineAny(start, end);
}


Foam::searchableSurface::volumeType
 Foam::triSurfaceMesh::getVolumeType
(
    const point& pt
) const
{
    // - use cached volume type per each tree node
    // - cheat conversion since same values
    return static_cast<volumeType>(tree().getVolumeType(pt));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //

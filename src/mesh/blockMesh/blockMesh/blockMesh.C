/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "blockMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::blockMesh::blockMesh::verboseOutput(false);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMesh::blockMesh(IOdictionary& dict)
:
    blockPointField_(dict.lookup("vertices")),
    scaleFactor_(1.0),
    topologyPtr_(createTopology(dict))
{
    calcMergeInfo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockMesh::~blockMesh()
{
    delete topologyPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockMesh::verbose(const bool on)
{
    verboseOutput = on;
}


const Foam::pointField& Foam::blockMesh::blockPointField() const
{
    return blockPointField_;
}


const Foam::polyMesh& Foam::blockMesh::topology() const
{
    if (!topologyPtr_)
    {
        FatalErrorIn("blockMesh::topology() const")
            << "topologyPtr_ not allocated"
            << exit(FatalError);
    }

    return *topologyPtr_;
}


Foam::scalar Foam::blockMesh::scaleFactor() const
{
    return scaleFactor_;
}


const Foam::pointField& Foam::blockMesh::points() const
{
    if (points_.empty())
    {
        createPoints();
    }

    return points_;
}


const Foam::cellShapeList& Foam::blockMesh::cells() const
{
    if (cells_.empty())
    {
        createCells();
    }

    return cells_;
}


const Foam::faceListList& Foam::blockMesh::patches() const
{
    if (patches_.empty())
    {
        createPatches();
    }

    return patches_;
}


Foam::wordList Foam::blockMesh::patchNames() const
{
    return topology().boundaryMesh().names();
}


Foam::wordList Foam::blockMesh::patchTypes() const
{
    return topology().boundaryMesh().types();
}


Foam::wordList Foam::blockMesh::patchPhysicalTypes() const
{
    return topology().boundaryMesh().physicalTypes();
}


Foam::label Foam::blockMesh::numZonedBlocks() const
{
    label num = 0;

    forAll(*this, blockI)
    {
        if (operator[](blockI).zoneName().size())
        {
            num++;
        }
    }

    return num;
}


void Foam::blockMesh::writeTopology(Ostream& os) const
{
    const pointField& pts = topology().points();

    forAll(pts, pI)
    {
        const point& pt = pts[pI];

        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }

    const edgeList& edges = topology().edges();

    forAll(edges, eI)
    {
        const edge& e = edges[eI];

        os << "l " << e.start() + 1 << ' ' << e.end() + 1 << endl;
    }
}

// ************************************************************************* //

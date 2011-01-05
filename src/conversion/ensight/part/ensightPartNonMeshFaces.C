/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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

#include "ensightPartNonMeshFaces.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
   defineTypeNameAndDebug(ensightPartNonMeshFaces, 0);
   addToRunTimeSelectionTable(ensightPart, ensightPartNonMeshFaces, istream);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::ensightPart::localPoints
Foam::ensightPartNonMeshFaces::calcLocalPoints() const
{
    localPoints ptList;
    ptList.list = identity(points_.size());
    ptList.nPoints = points_.size();
    return ptList;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightPartNonMeshFaces::ensightPartNonMeshFaces
(
    label partNumber,
    const string& partDescription,
    const faceList& faces,
    const pointField& points
)
:
    ensightPartFaces(partNumber, partDescription),
    faces_(faces),
    points_(points)
{
    binShapes(faces);
}


//- Construct as copy
Foam::ensightPartNonMeshFaces::ensightPartNonMeshFaces
(
    const ensightPartNonMeshFaces& part
)
:
    ensightPartFaces(part),
    faces_(part.faces_),
    points_(part.points_)
{}


//- Construct from Istream
Foam::ensightPartNonMeshFaces::ensightPartNonMeshFaces(Istream& is)
:
    ensightPartFaces(is),
    faces_(is),
    points_(is)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightPartNonMeshFaces::~ensightPartNonMeshFaces()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::ensightPartNonMeshFaces::writeConnectivity
(
    ensightGeoFile& os,
    const word& key,
    const labelList& idList,
    const labelList& pointMap
) const
{
    ensightPartFaces::writeConnectivity
    (
        os,
        key,
        faces_,
        idList,
        pointMap
    );
}


void Foam::ensightPartNonMeshFaces::writeGeometry(ensightGeoFile& os) const
{
    ensightPart::writeGeometry(os, points_);
}


// ************************************************************************* //

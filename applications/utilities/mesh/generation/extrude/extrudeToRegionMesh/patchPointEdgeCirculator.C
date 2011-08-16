/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "patchPointEdgeCirculator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::patchPointEdgeCirculator
Foam::patchPointEdgeCirculator::endConstIter
(
    *reinterpret_cast<primitiveFacePatch*>(0),  // primitiveFacePatch
    *reinterpret_cast<PackedBoolList*>(0),      // PackedBoolList
    -1,                                         // edgeID
    -1,                                         // index
    -1                                          // pointID
);


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<patchPointEdgeCirculator>& ip
)
{
    const patchPointEdgeCirculator& c = ip.t_;

    const pointField& pts = c.patch_.localPoints();
    const edge& e = c.patch_.edges()[c.edgeID_];
    label faceI = c.faceID();

    os  << "around point:" << c.pointID_
        << " coord:" << pts[c.pointID_]
        << " at edge:" << c.edgeID_
        << " coord:" << pts[e.otherVertex(c.pointID_)]
        << " in direction of face:" << faceI;

    if (faceI != -1)
    {
        os  << " fc:" << c.patch_.localFaces()[faceI].centre(pts);
    }
    return os;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "autoSnapDriver.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList>
Foam::labelList Foam::autoSnapDriver::getFacePoints
(
    const indirectPrimitivePatch& pp,
    const FaceList& faces
)
{
    // Could use PrimitivePatch & localFaces to extract points but might just
    // as well do it ourselves.

    boolList pointOnZone(pp.nPoints(), false);

    forAll(faces, i)
    {
        const face& f = faces[i];

        forAll(f, fp)
        {
            label meshPointI = f[fp];

            Map<label>::const_iterator iter =
                pp.meshPointMap().find(meshPointI);

            if (iter != pp.meshPointMap().end())
            {
                label pointI = iter();
                pointOnZone[pointI] = true;
            }
        }
    }

    return findIndices(pointOnZone, true);
}


// ************************************************************************* //

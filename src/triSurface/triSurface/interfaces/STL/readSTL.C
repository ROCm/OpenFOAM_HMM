/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "STLReader.H"
#include "mergePoints.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::triSurface::readSTL(const fileName& STLfileName, bool forceBinary)
{
    // Read in the values
    fileFormats::STLReader reader
    (
        STLfileName,
        (
            forceBinary
          ? fileFormats::STLCore::BINARY
          : fileFormats::STLCore::UNKNOWN
        )
    );

    // Stitch points
    labelList pointMap;
    label nUniquePoints = mergePoints
    (
        reader.points(),
        (
            // With the merge distance depending on the input format
            (reader.stlFormat() == fileFormats::STLCore::BINARY ? 10 : 100)
          * SMALL
        ),
        false,                  // verbose
        pointMap                // old to new point map
    );

    const pointField& readpts = reader.points();
    const labelList&  zoneIds = reader.zoneIds();

    pointField& pointLst = storedPoints();
    List<Face>& faceLst  = storedFaces();

    // Sizing
    pointLst.setSize(nUniquePoints);
    faceLst.setSize(zoneIds.size());

    // Assign points
    forAll(readpts, pointi)
    {
        pointLst[pointMap[pointi]] = readpts[pointi];
    }

    // Assign triangles
    label pointi = 0;
    forAll(faceLst, i)
    {
        Face& f = faceLst[i];

        f[0] = pointMap[pointi++];
        f[1] = pointMap[pointi++];
        f[2] = pointMap[pointi++];
        f.region() = zoneIds[i];
    }

    // Set patch names (and sizes)
    // - there is likely a more efficient means of doing this
    if (reader.stlFormat() == fileFormats::STLCore::ASCII)
    {
        const List<word>& names = reader.names();

        patches_.setSize(names.size());
        forAll(names, namei)
        {
            patches_[namei].name() = names[namei];
        }
        setDefaultPatches();
    }

    return true;
}

// ************************************************************************* //

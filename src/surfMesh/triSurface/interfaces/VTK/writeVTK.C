/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "triSurface.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// File-scope constant.
//
// TODO: make this run-time selectable (ASCII | BINARY)
// - Legacy mode only

static const Foam::foamVtkOutput::formatType fmtType =
    Foam::foamVtkOutput::formatType::LEGACY_ASCII;
    // Foam::foamVtkOutput::formatType::LEGACY_BASE64;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurface::writeVTK
(
    const bool writeSorted,
    std::ostream& os
) const
{
    autoPtr<foamVtkOutput::formatter> format =
        foamVtkOutput::newFormatter(os, fmtType);

    // Header
    foamVtkOutput::legacy::fileHeader
    (
        format(),
        "triSurface",
        vtkFileTag::POLY_DATA
    );

    const pointField& pts = this->points();
    const auto& faceLst = this->surfFaces();

    const label nFaces = faceLst.size();

    // Points
    foamVtkOutput::legacy::beginPoints(os, pts.size());

    foamVtkOutput::writeList(format(), pts);
    format().flush();

    // connectivity count (without additional storage)
    // is simply 3 * nFaces

    foamVtkOutput::legacy::beginPolys(os, nFaces, 3*nFaces);

    labelList faceMap;
    surfacePatchList patches(calcPatches(faceMap));

    const bool useFaceMap = (writeSorted && patches.size() > 1);

    if (useFaceMap)
    {
        label faceIndex = 0;
        for (const surfacePatch& patch : patches)
        {
            forAll(patch, i)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];

                format().write(3);  // The size prefix
                foamVtkOutput::writeList(format(), f);
            }
        }
        format().flush();


        // Write regions (zones) as CellData
        if (patches.size() > 1)
        {
            foamVtkOutput::legacy::dataHeader
            (
                os,
                vtkFileTag::CELL_DATA,
                nFaces,
                1  // Only one field
            );

            foamVtkOutput::legacy::intField
            (
                os,
                "region",
                1, // nComponent
                nFaces
            );

            faceIndex = 0;
            for (const surfacePatch& patch : patches)
            {
                forAll(patch, i)
                {
                    const Face& f = faceLst[faceMap[faceIndex++]];
                    format().write(f.region());
                }
            }
            format().flush();
        }
    }
    else
    {
        // No faceMap (unsorted)

        for (const Face& f : faceLst)
        {
            format().write(3);  // The size prefix
            foamVtkOutput::writeList(format(), f);
        }
        format().flush();


        // Write regions (zones) as CellData

        foamVtkOutput::legacy::dataHeader
        (
            os,
            vtkFileTag::CELL_DATA,
            faceLst.size(),
            1  // Only one field
        );

        foamVtkOutput::legacy::intField
        (
            os,
            "region",
            1, // nComponent
            faceLst.size()
        );

        for (const Face& f : faceLst)
        {
            format().write(f.region());
        }
        format().flush();
    }
}


// ************************************************************************* //

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
#include "STLCore.H"
#include "primitivePatch.H"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// A file-scope helper class to expose static member(s)
// This is a temporary measure and is expected to disappear in the future
struct triSurfaceSTLCore
:
    public Foam::fileFormats::STLCore
{
    using Foam::fileFormats::STLCore::writeBinaryHeader;
};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurface::writeSTLASCII(const bool writeSorted, Ostream& os) const
{
    labelList faceMap;

    surfacePatchList patches(calcPatches(faceMap));

    if (writeSorted)
    {
        label faceIndex = 0;
        forAll(patches, patchi)
        {
            // Print all faces belonging to this region
            const surfacePatch& patch = patches[patchi];

            os  << "solid " << patch.name() << endl;

            for
            (
                label patchFacei = 0;
                patchFacei < patch.size();
                patchFacei++
            )
            {
                const label facei = faceMap[faceIndex++];
                const labelledTri& f = (*this)[facei];

                // Write ASCII
                STLtriangle::write
                (
                    os,
                    faceNormals()[facei],
                    points()[f[0]],
                    points()[f[1]],
                    points()[f[2]]
                );
            }

            os  << "endsolid " << patch.name() << endl;
        }
    }
    else
    {
        // Get patch (=compact region) per face
        labelList patchIDs(size());
        forAll(patches, patchi)
        {
            label facei = patches[patchi].start();

            forAll(patches[patchi], i)
            {
                patchIDs[faceMap[facei++]] = patchi;
            }
        }

        label currentPatchi = -1;
        forAll(*this, facei)
        {
            if (currentPatchi != patchIDs[facei])
            {
                if (currentPatchi != -1)
                {
                    // Close previous solid
                    os  << "endsolid " << patches[currentPatchi].name() << nl;
                }
                currentPatchi = patchIDs[facei];
                os  << "solid " << patches[currentPatchi].name() << nl;
            }

            const labelledTri& f = (*this)[facei];

            // Write ASCII
            STLtriangle::write
            (
                os,
                faceNormals()[facei],
                points()[f[0]],
                points()[f[1]],
                points()[f[2]]
            );
        }

        if (currentPatchi != -1)
        {
            os  << "endsolid " << patches[currentPatchi].name() << nl;
        }
    }
}


void Foam::triSurface::writeSTLBINARY(std::ostream& os) const
{
    // Write the STL header
    triSurfaceSTLCore::writeBinaryHeader(os, this->size());

    forAll(*this, facei)
    {
        const labelledTri& f = (*this)[facei];

        // Write BINARY
        STLtriangle
        (
            faceNormals()[facei],
            points()[f[0]],
            points()[f[1]],
            points()[f[2]],
            f.region()
        ).write(os);
    }
}


// ************************************************************************* //

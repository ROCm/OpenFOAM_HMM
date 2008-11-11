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

#include "meshedSurface.H"
#include "mergePoints.H"
#include "triFace.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Remove badly degenerate faces, double faces.
void Foam::meshedSurface::cleanup(const bool verbose)
{
    // merge points (already done for STL, TRI)
    stitchFaces(SMALL, verbose);

    checkFaces(verbose);
    checkEdges(verbose);
}


bool Foam::meshedSurface::stitchFaces(const scalar tol, const bool verbose)
{
    pointField& pointLst = points();

    // Merge points
    labelList pointMap(pointLst.size());
    pointField newPoints(pointLst.size());

    bool hasMerged = mergePoints(pointLst, tol, verbose, pointMap, newPoints);

    if (!hasMerged)
    {
        return false;
    }

    if (verbose)
    {
        Info<< "meshedSurface::stitchFaces : Renumbering all faces"
            << endl;
    }

    // Set the coordinates to the merged ones
    pointLst.transfer(newPoints);

    List<FaceType>& faceLst = faces();

    // ensure we have at some patches, and they cover all the faces
    checkPatches();

    // Reset the point labels to the unique points array
    label oldFaceI = 0;
    label newFaceI = 0;
    forAll (patches_, patchI)
    {
        surfGroup& p = patches_[patchI];

        // adjust patch start
        p.start() = newFaceI;

        label patchEnd = oldFaceI + p.size();
        for (; oldFaceI < patchEnd; ++oldFaceI)
        {
            FaceType& f = faceLst[oldFaceI];
            forAll (f, fp)
            {
                f[fp] = pointMap[f[fp]];
            }

            if (f.collapse() >= 3)
            {
                if (newFaceI != oldFaceI)
                {
                    faceLst[newFaceI] = f;
                }
                newFaceI++;
            }
            else if (verbose)
            {
                Pout<< "meshedSurface::stitchFaces : "
                    << "Removing collapsed face " << oldFaceI << endl
                    << "    vertices   :" << f << endl;
            }
        }

        // adjust patch size
        p.size() = newFaceI - p.size();
    }

    if (newFaceI != faceLst.size())
    {
        if (verbose)
        {
            Pout<< "meshedSurface::stitchFaces : "
                << "Removed " << faceLst.size() - newFaceI
                << " faces" << endl;
        }
        faceLst.setSize(newFaceI);
    }


    // Merging points might have changed geometric factors
    MeshStorage::clearOut();

    return true;
}



// Remove badly degenerate faces and double faces.
void Foam::meshedSurface::checkFaces(const bool verbose)
{
    // Simple check on indices ok.
    const label maxPointI = points().size() - 1;

    List<FaceType>& faceLst = faces();

    // Phase 0: detect badly labelled faces
    forAll (faceLst, faceI)
    {
        const FaceType& f = faceLst[faceI];

        forAll (f, fp)
        {
            if (f[fp] < 0 || f[fp] > maxPointI)
            {
                FatalErrorIn("meshedSurface::checkFaces(bool)")
                    << "face " << f
                    << " uses point indices outside point range 0.."
                    << maxPointI
                    << exit(FatalError);
            }
        }
    }

    // ensure we have patches, and they cover all the faces
    checkPatches();

    // Phase 1: find and skip over invalid faces
    // Phase 2: pack
    const labelListList& fFaces = faceFaces();

    label oldFaceI = 0;
    label newFaceI = 0;
    forAll (patches_, patchI)
    {
        surfGroup& p = patches_[patchI];

        // correct the patch start
        p.start() = newFaceI;

        label patchEnd = oldFaceI + p.size();
        for (; oldFaceI < patchEnd; ++oldFaceI)
        {
            FaceType& f = faceLst[oldFaceI];

            // 'degenerate' face check
            if (f.collapse() >= 3)
            {
                // duplicate face check
                bool okay = true;
                const labelList& neighbours = fFaces[oldFaceI];

                // Check if faceNeighbours use same points as this face.
                // Note: discards normal information - sides of baffle are merged.
                forAll (neighbours, neighI)
                {
                    if (neighbours[neighI] <= oldFaceI)
                    {
                        // lower numbered faces already checked
                        continue;
                    }

                    const face& nei = faceLst[neighbours[neighI]];

                    if (f == nei)
                    {
                        okay = false;

                        if (verbose)
                        {
                            WarningIn
                            (
                                "meshedSurface::checkFaces(bool verbose)"
                            )   << "faces share the same vertices:\n"
                                << "    face 1 :" << oldFaceI << endl;
                            // printFace(Warning, "    ", f, points());

                            Warning
                                << endl
                                << "    face 2 :"
                                << neighbours[neighI] << endl;
                            // printFace(Warning, "    ", nei, points());
                        }

                        break;
                    }
                }

                if (okay)
                {
                    if (newFaceI != oldFaceI)
                    {
                        faceLst[newFaceI] = f;
                    }
                    newFaceI++;
                }
            }
            else if (verbose)
            {
                WarningIn
                (
                    "meshedSurface::checkFaces(bool verbose)"
                )   << "face " << oldFaceI
                    << " has fewer than three unique vertices:\n";
                // printTriangle(Warning, "    ", f, points());
            }
        }

        // adjust patch size
        p.size() = newFaceI - p.start();
    }

    if (newFaceI < faceLst.size())
    {
        if (verbose)
        {
            WarningIn
            (
                "meshedSurface::checkFaces(bool verbose)"
            )   << "Removed " << faceLst.size() - newFaceI
                << " illegal faces." << endl;
        }
        faceLst.setSize(newFaceI);

        // Topology can change because of renumbering
        MeshStorage::clearOut();
    }
}


Foam::label Foam::meshedSurface::triangulate()
{
    label nTri = 0;
    List<FaceType>& faceLst = faces();

    // determine how many triangles are needed
    forAll (faceLst, faceI)
    {
        nTri += faceLst[faceI].size() - 2;
    }

    // nothing to do
    if (nTri <= faceLst.size())
    {
        return 0;
    }

    List<FaceType> newFaces(nTri);

    // note the number of *additional* faces
    nTri -= faceLst.size();

    // Reset the point labels to the unique points array
    label oldFaceI = 0;
    label newFaceI = 0;
    forAll (patches_, patchI)
    {
        surfGroup& p = patches_[patchI];

        // adjust patch start
        p.start() = newFaceI;

        label patchEnd = oldFaceI + p.size();
        for (; oldFaceI < patchEnd; ++oldFaceI)
        {
            const FaceType& f = faceLst[oldFaceI];
            triFace fTri;

            // Do simple face triangulation around f[0].
            // we could also use face::triangulation
            fTri[0] = f[0];
            for (label fp = 1; fp < f.size() - 1; ++fp)
            {
                label fp1 = (fp + 1) % f.size();

                fTri[1] = f[fp];
                fTri[2] = f[fp1];

                newFaces[newFaceI++] = fTri;
            }
        }

        // adjust patch size
        p.size() = newFaceI - p.start();
    }

    faceLst.transfer(newFaces);

    return nTri;
}

// ************************************************************************* //

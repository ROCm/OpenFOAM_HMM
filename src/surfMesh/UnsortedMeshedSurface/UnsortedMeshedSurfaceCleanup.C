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

#include "UnsortedMeshedSurface.H"
#include "mergePoints.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Remove badly degenerate faces, double faces.
template<class Face>
void Foam::UnsortedMeshedSurface<Face>::cleanup(const bool verbose)
{
    // merge points (already done for STL, TRI)
    stitchFaces(SMALL, verbose);

    checkFaces(verbose);
    ParentType::checkEdges(verbose);
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::stitchFaces
(
    const scalar tol,
    const bool verbose
)
{
    pointField& pointLst = this->storedPoints();

    // Merge points
    labelList  pointMap(pointLst.size());
    pointField newPoints(pointLst.size());

    bool hasMerged = mergePoints(pointLst, tol, verbose, pointMap, newPoints);

    if (!hasMerged)
    {
        return false;
    }

    if (verbose)
    {
        Info<< "UnsortedMeshedSurface::stitchFaces : Renumbering all faces"
            << endl;
    }

    // Set the coordinates to the merged ones
    pointLst.transfer(newPoints);

    List<Face>& faceLst = this->storedFaces();

    // Reset the point labels to the unique points array
    label newFaceI = 0;
    forAll(faceLst, faceI)
    {
        Face& f = faceLst[faceI];
        forAll(f, fp)
        {
            f[fp] = pointMap[f[fp]];
        }

        if (f.collapse() >= 3)
        {
            if (newFaceI != faceI)
            {
                faceLst[newFaceI] = f;
                regions_[newFaceI] = regions_[faceI];
            }
            newFaceI++;
        }
        else if (verbose)
        {
            Pout<< "UnsortedMeshedSurface::stitchFaces : "
                << "Removing collapsed face " << faceI << endl
                << "    vertices   :" << f << endl;
        }
    }

    if (newFaceI != faceLst.size())
    {
        if (verbose)
        {
            Pout<< "UnsortedMeshedSurface::stitchFaces : "
                << "Removed " << faceLst.size() - newFaceI
                << " faces" << endl;
        }
        faceLst.setSize(newFaceI);
        regions_.setSize(newFaceI);
    }

    // Merging points might have changed geometric factors
    ParentType::clearOut();

    return true;
}


// Remove badly degenerate faces and double faces.
template<class Face>
void Foam::UnsortedMeshedSurface<Face>::checkFaces(const bool verbose)
{
    // Simple check on indices ok.
    const label maxPointI = this->points().size() - 1;

    List<Face>& faceLst = this->storedFaces();

    // Phase 0: detect badly labelled faces
    forAll(faceLst, faceI)
    {
        const Face& f = faceLst[faceI];

        forAll(f, fp)
        {
            if (f[fp] < 0 || f[fp] > maxPointI)
            {
                FatalErrorIn("UnsortedMeshedSurface::checkFaces(bool)")
                    << "face " << f
                    << " uses point indices outside point range 0.."
                    << maxPointI
                    << exit(FatalError);
            }
        }
    }

    // Phase 1: mark invalid faces
    // Phase 1: pack
    // Done to keep numbering constant in phase 1
    const labelListList& fFaces = ParentType::faceFaces();
    label newFaceI = 0;

    forAll(faceLst, faceI)
    {
        Face& f = faceLst[faceI];

        // avoid degenerate faces
        if (f.collapse() >= 3)
        {
            // duplicate face check
            bool okay = true;
            const labelList& neighbours = fFaces[faceI];

            // Check if faceNeighbours use same points as this face.
            // Note: discards normal information - sides of baffle are merged.
            forAll(neighbours, neighI)
            {
                if (neighbours[neighI] <= faceI)
                {
                    // lower numbered faces already checked
                    continue;
                }

                const Face& nei = faceLst[neighbours[neighI]];

                if (f == nei)
                {
                    okay = false;

                    if (verbose)
                    {
                        WarningIn
                        (
                            "UnsortedMeshedSurface::checkFaces(bool verbose)"
                        )   << "faces share the same vertices:\n"
                            << "    face 1 :" << faceI << endl;
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
                if (newFaceI != faceI)
                {
                    faceLst[newFaceI] = f;
                    regions_[newFaceI] = regions_[faceI];
                }
                newFaceI++;
            }
        }
        else if (verbose)
        {
            WarningIn
            (
                "UnsortedMeshedSurface::checkFaces(bool verbose)"
            )   << "face " << faceI
                << " does not at least three unique vertices:\n";
            // printFace(Warning, "    ", f, points());
        }
    }

    if (newFaceI < faceLst.size())
    {
        if (verbose)
        {
            WarningIn
            (
                "UnsortedMeshedSurface::checkFaces(bool verbose)"
            )   << "Removed " << faceLst.size() - newFaceI
                << " illegal faces." << endl;
        }
        faceLst.setSize(newFaceI);
        regions_.setSize(newFaceI);

        // Topology can change because of renumbering
        ParentType::clearOut();
    }
}


template<class Face>
Foam::label Foam::UnsortedMeshedSurface<Face>::triangulate()
{
    label nTri = 0;
    List<Face>& faceLst = this->storedFaces();

    // determine how many triangles are needed
    forAll(faceLst, faceI)
    {
        nTri += faceLst[faceI].size() - 2;
    }

    // nothing to do
    if (nTri <= faceLst.size())
    {
        return 0;
    }

    List<Face>  newFaces(nTri);
    List<label> newRegions(nTri);

    // note the number of *additional* faces
    nTri -= faceLst.size();

    // Reset the point labels to the unique points array
    label newFaceI = 0;
    forAll(faceLst, faceI)
    {
        const Face& f = faceLst[faceI];
        triFace fTri;

        // Do simple face triangulation around f[0].
        // we could also use face::triangulation
        fTri[0] = f[0];
        for (label fp = 1; fp < f.size() - 1; ++fp)
        {
            label fp1 = (fp + 1) % f.size();

            fTri[1] = f[fp];
            fTri[2] = f[fp1];

            newFaces[newFaceI] = fTri;
            newRegions[newFaceI] = regions_[faceI];
            newFaceI++;
        }
    }

    faceLst.transfer(newFaces);
    regions_.transfer(newRegions);

    return nTri;
}

// ************************************************************************* //

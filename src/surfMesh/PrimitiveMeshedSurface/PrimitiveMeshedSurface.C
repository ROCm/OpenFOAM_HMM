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

#include "PrimitiveMeshedSurface.H"
#include "mergePoints.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
inline bool Foam::PrimitiveMeshedSurface<Face>::isTri()
{
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::PrimitiveMeshedSurface<Face>::PrimitiveMeshedSurface()
:
    ParentType(List<Face>(), pointField())
{}


template<class Face>
Foam::PrimitiveMeshedSurface<Face>::PrimitiveMeshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst
)
:
    ParentType(List<Face>(), pointField())
{
    reset(pointLst, faceLst);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::PrimitiveMeshedSurface<Face>::~PrimitiveMeshedSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::PrimitiveMeshedSurface<Face>::clear()
{
    ParentType::clearOut();

    storedPoints().clear();
    storedFaces().clear();
}


template<class Face>
void Foam::PrimitiveMeshedSurface<Face>::movePoints(const pointField& newPoints)
{
    // Remove all geometry dependent data
    ParentType::clearTopology();

    // Adapt for new point position
    ParentType::movePoints(newPoints);

    // Copy new points
    storedPoints() = newPoints;
}


template<class Face>
void Foam::PrimitiveMeshedSurface<Face>::scalePoints(const scalar& scaleFactor)
{
    // avoid bad scaling
    if (scaleFactor > 0 && scaleFactor != 1.0)
    {
        // Remove all geometry dependent data
        ParentType::clearTopology();

        // Adapt for new point position
        ParentType::movePoints(pointField());

        storedPoints() *= scaleFactor;
    }
}


template<class Face>
void Foam::PrimitiveMeshedSurface<Face>::reset
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst
)
{
    ParentType::clearOut();

    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (&pointLst)
    {
        storedPoints().transfer(pointLst());
    }

    if (&faceLst)
    {
        storedFaces().transfer(faceLst());
    }
}


// Remove badly degenerate faces, double faces.
template<class Face>
void Foam::PrimitiveMeshedSurface<Face>::cleanup(const bool verbose)
{
    // merge points (already done for STL, TRI)
    stitchFaces(SMALL, verbose);

    checkFaces(verbose);
    this->checkEdges(verbose);
}


template<class Face>
bool Foam::PrimitiveMeshedSurface<Face>::stitchFaces
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
        Info<< "PrimitiveMeshedSurface::stitchFaces : Renumbering all faces"
            << endl;
    }

    // Set the coordinates to the merged ones
    pointLst.transfer(newPoints);

    List<Face>& faceLst = this->storedFaces();

    List<label> faceMap(faceLst.size());

    // Reset the point labels to the unique points array
    label newFaceI = 0;
    forAll(faceLst, faceI)
    {
        Face& f = faceLst[faceI];
        forAll(f, fp)
        {
            f[fp] = pointMap[f[fp]];
        }

        // for extra safety: collapse face as well
        if (f.collapse() >= 3)
        {
            if (newFaceI != faceI)
            {
                faceLst[newFaceI] = f;
            }
            faceMap[newFaceI] = faceI;
            newFaceI++;
        }
        else if (verbose)
        {
            Pout<< "PrimitiveMeshedSurface::stitchFaces : "
                << "Removing collapsed face " << faceI << endl
                << "    vertices   :" << f << endl;
        }
    }
    pointMap.clear();

    if (newFaceI != faceLst.size())
    {
        if (verbose)
        {
            Pout<< "PrimitiveMeshedSurface::stitchFaces : "
                << "Removed " << faceLst.size() - newFaceI
                << " faces" << endl;
        }
        faceLst.setSize(newFaceI);
        remapRegions(faceMap);
    }
    faceMap.clear();

    // Merging points might have changed geometric factors
    ParentType::clearOut();
    return true;
}


// Remove badly degenerate faces and double faces.
template<class Face>
bool Foam::PrimitiveMeshedSurface<Face>::checkFaces
(
    const bool verbose
)
{
    bool changed = false;
    List<Face>& faceLst = this->storedFaces();

    List<label> faceMap(faceLst.size());

    label newFaceI = 0;
    // Detect badly labelled faces and mark degenerate faces
    const label maxPointI = this->points().size() - 1;
    forAll(faceLst, faceI)
    {
        Face& f = faceLst[faceI];

        // avoid degenerate faces
        if (f.collapse() >= 3)
        {
            forAll(f, fp)
            {
                if (f[fp] < 0 || f[fp] > maxPointI)
                {
                    FatalErrorIn("PrimitiveMeshedSurface::checkFaces(bool)")
                        << "face " << f
                        << " uses point indices outside point range 0.."
                    << maxPointI
                        << exit(FatalError);
                }
            }

            faceMap[faceI] = faceI;
            newFaceI++;
        }
        else
        {
            // mark as bad face
            faceMap[faceI] = -1;

            changed = true;
            if (verbose)
            {
                WarningIn
                (
                    "PrimitiveMeshedSurface::checkFaces(bool verbose)"
                )   << "face[" << faceI << "] = " << f
                    << " does not have three unique vertices" << endl;
            }
        }
    }

    // Detect doubled faces
    // do not touch the faces
    const labelListList& fFaces = this->faceFaces();
    newFaceI = 0;
    forAll(faceLst, faceI)
    {
        // skip already collapsed faces:
        if (faceMap[faceI] < 0)
        {
            continue;
        }

        const Face& f = faceLst[faceI];

        // duplicate face check
        bool okay = true;
        const labelList& neighbours = fFaces[faceI];

        // Check if faceNeighbours use same points as this face.
        // Note: discards normal information - sides of baffle are merged.
        forAll(neighbours, neighI)
        {
            const label neiFaceI = neighbours[neighI];

            if (neiFaceI <= faceI || faceMap[neiFaceI] < 0)
            {
                // lower numbered faces already checked
                // skip neighbours that are themselves collapsed
                continue;
            }

            const Face& nei = faceLst[neiFaceI];

            if (f == nei)
            {
                okay = false;

                if (verbose)
                {
                    WarningIn
                    (
                        "PrimitiveMeshedSurface::checkFaces(bool verbose)"
                    )   << "faces share the same vertices:" << nl
                        << "    face[" << faceI << "] : " << f << nl
                        << "    face[" << neiFaceI << "] : " << nei << endl;
                    // printFace(Warning, "    ", f, points());
                    // printFace(Warning, "    ", nei, points());
                }

                break;
            }
        }

        if (okay)
        {
            faceMap[faceI] = faceI;
            newFaceI++;
        }
        else
        {
            faceMap[faceI] = -1;
        }
    }

    // Phase 1: pack
    // Done to keep numbering constant in phase 1

    if (changed || newFaceI < faceLst.size())
    {
        changed = true;

        if (verbose)
        {
            WarningIn
            (
                "PrimitiveMeshedSurface::checkFaces(bool verbose)"
            )   << "Removed " << faceLst.size() - newFaceI
                << " illegal faces." << endl;
        }

        // compress the face list
        newFaceI = 0;
        forAll(faceLst, faceI)
        {
            if (faceMap[faceI] >= 0)
            {
                if (newFaceI != faceI)
                {
                    faceLst[newFaceI] = faceLst[faceI];
                }
                faceMap[newFaceI] = faceI;
                newFaceI++;
            }
        }

        faceLst.setSize(newFaceI);
        remapRegions(faceMap);
    }
    faceMap.clear();

    // Topology can change because of renumbering
    ParentType::clearOut();
    return changed;
}


template<class Face>
Foam::label Foam::PrimitiveMeshedSurface<Face>::triangulate()
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
    List<label> faceMap(nTri);

    // remember the number of *additional* faces
    nTri -= faceLst.size();

    label newFaceI = 0;
    forAll(faceLst, faceI)
    {
        const Face& f = faceLst[faceI];
        triFace fTri;

        // Do simple face triangulation around f[0].
        // we could also use face::triangulation, but that requires points
        // and doesn't currently template nicely
        fTri[0] = f[0];
        for (label fp = 1; fp < f.size() - 1; ++fp)
        {
            label fp1 = (fp + 1) % f.size();

            fTri[1] = f[fp];
            fTri[2] = f[fp1];

            newFaces[newFaceI] = fTri;
            faceMap[newFaceI] = faceI;
            newFaceI++;
        }
    }

    faceLst.transfer(newFaces);
    remapRegions(faceMap);
    faceMap.clear();

    // Topology can change because of renumbering
    ParentType::clearOut();
    return nTri;
}


// dummy implementation to avoid a pure virtual class
template<class Face>
void Foam::PrimitiveMeshedSurface<Face>::remapRegions(List<label>& faceMap)
{
    faceMap.clear();
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

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

#include "BasicMeshedSurface.H"
#include "boundBox.H"
#include "mergePoints.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
inline bool Foam::BasicMeshedSurface<Face>::isTri()
{
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::BasicMeshedSurface<Face>::BasicMeshedSurface()
:
    ParentType(List<Face>(), pointField())
{}


template<class Face>
Foam::BasicMeshedSurface<Face>::BasicMeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst
)
:
    ParentType(List<Face>(), pointField())
{
    reset(pointLst, faceLst);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::BasicMeshedSurface<Face>::~BasicMeshedSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::BasicMeshedSurface<Face>::clear()
{
    ParentType::clearOut();

    storedPoints().clear();
    storedFaces().clear();
}


template<class Face>
void Foam::BasicMeshedSurface<Face>::movePoints(const pointField& newPoints)
{
    // Remove all geometry dependent data
    ParentType::clearTopology();

    // Adapt for new point position
    ParentType::movePoints(newPoints);

    // Copy new points
    storedPoints() = newPoints;
}


template<class Face>
void Foam::BasicMeshedSurface<Face>::scalePoints(const scalar& scaleFactor)
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
void Foam::BasicMeshedSurface<Face>::reset
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst
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


template<class Face>
void Foam::BasicMeshedSurface<Face>::reset
(
    const Xfer<List<point> >& pointLst,
    const Xfer<List<Face> >& faceLst
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
void Foam::BasicMeshedSurface<Face>::cleanup(const bool verbose)
{
    // merge points (already done for STL, TRI)
    stitchFaces(SMALL, verbose);

    checkFaces(verbose);
    this->checkTopology(verbose);
}


template<class Face>
bool Foam::BasicMeshedSurface<Face>::stitchFaces
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
        Info<< "BasicMeshedSurface::stitchFaces : Renumbering all faces"
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
            Pout<< "BasicMeshedSurface::stitchFaces : "
                << "Removing collapsed face " << faceI << endl
                << "    vertices   :" << f << endl;
        }
    }
    pointMap.clear();

    if (newFaceI != faceLst.size())
    {
        if (verbose)
        {
            Pout<< "BasicMeshedSurface::stitchFaces : "
                << "Removed " << faceLst.size() - newFaceI
                << " faces" << endl;
        }
        faceLst.setSize(newFaceI);
        remapFaces(faceMap);
    }
    faceMap.clear();

    // Merging points might have changed geometric factors
    ParentType::clearOut();
    return true;
}


// Remove badly degenerate faces and double faces.
template<class Face>
bool Foam::BasicMeshedSurface<Face>::checkFaces
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
                    FatalErrorIn("BasicMeshedSurface::checkFaces(bool)")
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
                    "BasicMeshedSurface::checkFaces(bool verbose)"
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
                        "BasicMeshedSurface::checkFaces(bool verbose)"
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
                "BasicMeshedSurface::checkFaces(bool verbose)"
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
        remapFaces(faceMap);
    }
    faceMap.clear();

    // Topology can change because of renumbering
    ParentType::clearOut();
    return changed;
}


template<class Face>
Foam::label Foam::BasicMeshedSurface<Face>::triangulate()
{
    return triangulate
    (
        const_cast<List<label>&>(List<label>::null())
    );
}


template<class Face>
Foam::label Foam::BasicMeshedSurface<Face>::triangulate
(
    List<label>& faceMapOut
)
{
    label nTri = 0;
    label maxTri = 0;  // the maximum number of triangles for any single face
    List<Face>& faceLst = this->storedFaces();

    // determine how many triangles will be needed
    forAll(faceLst, faceI)
    {
        const label n = faceLst[faceI].nTriangles();
        if (maxTri < n)
        {
            maxTri = n;
        }
        nTri += n;
    }

    // nothing to do
    if (nTri <= faceLst.size())
    {
        if (&faceMapOut)
        {
            faceMapOut.clear();
        }
        return 0;
    }

    List<Face>  newFaces(nTri);
    List<label> faceMap;

    // reuse storage from optional faceMap
    if (&faceMapOut)
    {
        faceMap.transfer(faceMapOut);
    }
    faceMap.setSize(nTri);

    // remember the number of *additional* faces
    nTri -= faceLst.size();

    if (this->points().empty())
    {
        // triangulate without points
        // simple face triangulation around f[0]
        label newFaceI = 0;
        forAll(faceLst, faceI)
        {
            const Face& f = faceLst[faceI];

            for (label fp = 1; fp < f.size() - 1; ++fp)
            {
                label fp1 = (fp + 1) % f.size();

                newFaces[newFaceI] = triFace(f[0], f[fp], f[fp1]);
                faceMap[newFaceI] = faceI;
                newFaceI++;
            }
        }
    }
    else
    {
        // triangulate with points
        List<face> tmpTri(maxTri);

        label newFaceI = 0;
        forAll(faceLst, faceI)
        {
            // 'face' not '<Face>'
            const face& f = faceLst[faceI];

            label nTmp;
            f.triangles(this->points(), nTmp, tmpTri);
            for (label triI = 0; triI < nTmp; triI++)
            {
                newFaces[newFaceI] = Face
                (
                    static_cast<UList<label>&>(tmpTri[triI])
                );
                faceMap[newFaceI] = faceI;
                newFaceI++;
            }
        }
    }

    faceLst.transfer(newFaces);
    remapFaces(faceMap);

    // optionally return the faceMap
    if (&faceMapOut)
    {
        faceMapOut.transfer(faceMap);
    }
    faceMap.clear();

    // Topology can change because of renumbering
    ParentType::clearOut();
    return nTri;
}


// dummy implementation to avoid a pure virtual class
template<class Face>
void Foam::BasicMeshedSurface<Face>::remapFaces(const UList<label>&)
{
}


template<class Face>
void Foam::BasicMeshedSurface<Face>::writeStats(Ostream& os) const
{
    os  << "points      : " << this->points().size() << nl;
    if (this->isTri())
    {
        os << "triangles   : " << this->size() << nl;
    }
    else
    {
        label nTri = 0;
        label nQuad = 0;
        forAll(*this, i)
        {
            const label n = this->operator[](i).size();

            if (n == 3)
            {
                nTri++;
            }
            else if (n == 4)
            {
                nQuad++;
            }
        }

        os  << "faces       : " << this->size()
            << "  (tri:" << nTri << " quad:" << nQuad 
            << " poly:" << (this->size() - nTri - nQuad ) << ")" << nl;
    }

    os  << "boundingBox : " << boundBox(this->points()) << endl;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

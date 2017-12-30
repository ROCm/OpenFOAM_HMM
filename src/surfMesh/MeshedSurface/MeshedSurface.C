/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "mergePoints.H"
#include "Time.H"
#include "ListOps.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "surfMesh.H"
#include "primitivePatch.H"
#include "faceTraits.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Face>
Foam::wordHashSet Foam::MeshedSurface<Face>::readTypes()
{
    return wordHashSet(*fileExtensionConstructorTablePtr_);
}

template<class Face>
Foam::wordHashSet Foam::MeshedSurface<Face>::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
bool Foam::MeshedSurface<Face>::canReadType
(
    const word& ext,
    const bool verbose
)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        readTypes() | FriendType::readTypes(),
        ext,
        verbose,
        "reading"
   );
}


template<class Face>
bool Foam::MeshedSurface<Face>::canWriteType
(
    const word& ext,
    const bool verbose
)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        writeTypes() | ProxyType::writeTypes(),
        ext,
        verbose,
        "writing"
    );
}


template<class Face>
bool Foam::MeshedSurface<Face>::canRead
(
    const fileName& name,
    const bool verbose
)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }
    return canReadType(ext, verbose);
}


template<class Face>
void Foam::MeshedSurface<Face>::write
(
    const fileName& name,
    const MeshedSurface<Face>& surf,
    const dictionary& options
)
{
    write(name, name.ext(), surf, options);
}


template<class Face>
void Foam::MeshedSurface<Face>::write
(
    const fileName& name,
    const word& ext,
    const MeshedSurface<Face>& surf,
    const dictionary& options
)
{
    if (debug)
    {
        InfoInFunction << "Writing to " << name << endl;
    }

    auto mfIter = writefileExtensionMemberFunctionTablePtr_->cfind(ext);

    if (!mfIter.found())
    {
        // No direct writer, delegate to proxy if possible
        const wordHashSet& delegate = ProxyType::writeTypes();

        if (delegate.found(ext))
        {
            MeshedSurfaceProxy<Face>(surf).write(name, ext, options);
        }
        else
        {
            FatalErrorInFunction
                << "Unknown file extension " << ext << nl << nl
                << "Valid types:" << nl
                << flatOutput((delegate | writeTypes()).sortedToc()) << nl
                << exit(FatalError);
        }
    }
    else
    {
        mfIter()(name, surf, options);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface()
:
    ParentType(List<Face>(), pointField()),
    zones_()
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face>>& faceLst,
    const Xfer<surfZoneList>& zoneLst
)
:
    MeshedSurface<Face>()
{
    reset(pointLst, faceLst, zoneLst);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face>>& faceLst,
    const labelUList& zoneSizes,
    const UList<word>& zoneNames
)
:
    MeshedSurface<Face>()
{
    reset(pointLst, faceLst, Xfer<surfZoneList>());

    if (zoneSizes.size())
    {
        if (zoneNames.size())
        {
            addZones(zoneSizes, zoneNames);
        }
        else
        {
            addZones(zoneSizes);
        }
    }
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const MeshedSurface<Face>& surf
)
:
    ParentType(surf.surfFaces(), surf.points()),
    zones_(surf.surfZones())
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const UnsortedMeshedSurface<Face>& surf
)
:
    ParentType(List<Face>(), surf.points())
{
    labelList faceMap;
    this->storedZones() = surf.sortedZones(faceMap);

    const List<Face>& origFaces = surf;
    List<Face> newFaces(origFaces.size());

    forAll(newFaces, facei)
    {
        newFaces[faceMap[facei]] = origFaces[facei];
    }

    this->storedFaces().transfer(newFaces);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const surfMesh& mesh)
:
    MeshedSurface<Face>()
{
    // same face type as surfMesh
    MeshedSurface<face> surf
    (
        xferCopy(mesh.points()),
        xferCopy(mesh.faces()),
        xferCopy(mesh.surfZones())
    );

    this->transcribe(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const polyBoundaryMesh& bMesh,
    const bool useGlobalPoints
)
:
    MeshedSurface<Face>()
{
    const polyMesh& mesh = bMesh.mesh();
    const polyPatchList& bPatches = bMesh;

    // Get a single patch for all boundaries
    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nFaces() - mesh.nInternalFaces(),
            mesh.nInternalFaces()
        ),
        mesh.points()
    );

    // use global/local points:
    const pointField& bPoints =
    (
        useGlobalPoints ? mesh.points() : allBoundary.localPoints()
    );

    // global/local face addressing:
    const List<Face>& bFaces =
    (
        useGlobalPoints ? allBoundary : allBoundary.localFaces()
    );


    // create zone list
    surfZoneList newZones(bPatches.size());

    label startFacei = 0;
    label nZone = 0;
    forAll(bPatches, patchi)
    {
        const polyPatch& p = bPatches[patchi];

        if (p.size())
        {
            newZones[nZone] = surfZone
            (
                p.name(),
                p.size(),
                startFacei,
                nZone
            );

            ++nZone;
            startFacei += p.size();
        }
    }

    newZones.setSize(nZone);

    // same face type as the polyBoundaryMesh
    MeshedSurface<face> surf
    (
        xferCopy(bPoints),
        xferCopy(bFaces),
        xferMove(newZones)
    );

    this->transcribe(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const fileName& name, const word& ext)
:
    MeshedSurface<Face>()
{
    read(name, ext);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const fileName& name)
:
    MeshedSurface<Face>()
{
    read(name);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(Istream& is)
:
    MeshedSurface<Face>()
{
    read(is);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Time& t,
    const word& surfName
)
:
    MeshedSurface<Face>()
{
    surfMesh mesh
    (
        IOobject
        (
            "dummyName",
            t.timeName(),
            t,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        ),
        surfName
    );

    // same face type as surfMesh
    MeshedSurface<face> surf
    (
        xferMove(mesh.storedPoints()),
        xferMove(mesh.storedFaces()),
        xferMove(mesh.storedZones())
    );

    this->transcribe(surf);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<MeshedSurface<Face>>& surf
)
:
    MeshedSurface<Face>()
{
    transfer(surf());
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<UnsortedMeshedSurface<Face>>& surf
)
:
    MeshedSurface<Face>()
{
    transfer(surf());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurface<Face>::~MeshedSurface()
{
    clear();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::remapFaces
(
    const labelUList& faceMap
)
{
    if (isNull(faceMap) || faceMap.empty())
    {
        return;
    }

    surfZoneList& zones = storedZones();

    if (zones.size() == 1)
    {
        zones[0].size() = faceMap.size();   // Single zone case is trivial
        return;
    }

    // Recalculate the zone start/size
    label newFacei = 0;
    label origEndI = 0;

    for (surfZone& zone : zones)
    {
        // Adjust zone start
        zone.start() = newFacei;
        origEndI += zone.size();

        for (label facei = newFacei; facei < faceMap.size(); ++facei)
        {
            if (faceMap[facei] < origEndI)
            {
                ++newFacei;
            }
            else
            {
                break;
            }
        }

        // Adjust zone size
        zone.size() = newFacei - zone.start();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::clear()
{
    ParentType::clearOut();

    storedPoints().clear();
    storedFaces().clear();
    storedZones().clear();
}


template<class Face>
void Foam::MeshedSurface<Face>::movePoints(const pointField& newPoints)
{
    // Changes areas, normals etc.
    ParentType::clearGeom();

    // Adapt for new point position
    ParentType::movePoints(newPoints);

    // Copy new points
    storedPoints() = newPoints;
}


template<class Face>
void Foam::MeshedSurface<Face>::scalePoints(const scalar scaleFactor)
{
    // Avoid bad scaling
    if (scaleFactor > 0 && scaleFactor != 1.0)
    {
        // Changes areas, normals etc.
        ParentType::clearGeom();

        pointField newPoints(scaleFactor*this->points());

        // Adapt for new point position
        ParentType::movePoints(newPoints);

        storedPoints() = newPoints;
    }
}


template<class Face>
void Foam::MeshedSurface<Face>::reset
(
    const Xfer<MeshedSurface<Face>>& surf
)
{
    transfer(surf());
}


template<class Face>
void Foam::MeshedSurface<Face>::reset
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face>>& faceLst,
    const Xfer<surfZoneList>& zoneLst
)
{
    ParentType::clearOut();

    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (notNull(pointLst))
    {
        storedPoints().transfer(pointLst());
    }

    if (notNull(faceLst))
    {
        storedFaces().transfer(faceLst());
    }

    if (notNull(zoneLst))
    {
        storedZones().transfer(zoneLst());
    }
}


template<class Face>
void Foam::MeshedSurface<Face>::reset
(
    const Xfer<List<point>>& pointLst,
    const Xfer<List<Face>>& faceLst,
    const Xfer<surfZoneList>& zoneLst
)
{
    ParentType::clearOut();

    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (notNull(pointLst))
    {
        storedPoints().transfer(pointLst());
    }

    if (notNull(faceLst))
    {
        storedFaces().transfer(faceLst());
    }

    if (notNull(zoneLst))
    {
        storedZones().transfer(zoneLst());
    }
}


// Remove badly degenerate faces, double faces.
template<class Face>
void Foam::MeshedSurface<Face>::cleanup(const bool verbose)
{
    // Merge points (already done for STL, TRI)
    stitchFaces(SMALL, verbose);

    checkFaces(verbose);
    this->checkTopology(verbose);
}


template<class Face>
bool Foam::MeshedSurface<Face>::stitchFaces
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
        InfoInFunction<< "Renumbering all faces" << endl;
    }

    // Set the coordinates to the merged ones
    pointLst.transfer(newPoints);

    List<Face>& faceLst = this->storedFaces();

    List<label> faceMap(faceLst.size());

    // Reset the point labels to the unique points array
    label newFacei = 0;
    forAll(faceLst, facei)
    {
        Face& f = faceLst[facei];
        forAll(f, fp)
        {
            f[fp] = pointMap[f[fp]];
        }

        // for extra safety: collapse face as well
        if (f.collapse() >= 3)
        {
            if (newFacei != facei)
            {
                faceLst[newFacei] = f;
            }
            faceMap[newFacei] = facei;
            ++newFacei;
        }
        else if (verbose)
        {
            Pout<< "MeshedSurface::stitchFaces : "
                << "Removing collapsed face " << facei << endl
                << "    vertices   :" << f << endl;
        }
    }
    pointMap.clear();

    if (newFacei != faceLst.size())
    {
        if (verbose)
        {
            Pout<< "MeshedSurface::stitchFaces : "
                << "Removed " << faceLst.size() - newFacei
                << " faces" << endl;
        }
        faceLst.setSize(newFacei);
        faceMap.setSize(newFacei);
        remapFaces(faceMap);
    }
    faceMap.clear();

    // Merging points might have changed geometric factors
    ParentType::clearOut();
    return true;
}


// Remove badly degenerate faces and double faces.
template<class Face>
bool Foam::MeshedSurface<Face>::checkFaces
(
    const bool verbose
)
{
    bool changed = false;
    List<Face>& faceLst = this->storedFaces();

    List<label> faceMap(faceLst.size());

    label newFacei = 0;
    // Detect badly labelled faces and mark degenerate faces
    const label maxPointi = this->points().size() - 1;
    forAll(faceLst, facei)
    {
        Face& f = faceLst[facei];

        // avoid degenerate faces
        if (f.collapse() >= 3)
        {
            forAll(f, fp)
            {
                if (f[fp] < 0 || f[fp] > maxPointi)
                {
                    FatalErrorInFunction
                        << "face " << f
                        << " uses point indices outside point range 0.."
                    << maxPointi
                        << exit(FatalError);
                }
            }

            faceMap[facei] = facei;
            ++newFacei;
        }
        else
        {
            // mark as bad face
            faceMap[facei] = -1;

            changed = true;
            if (verbose)
            {
                WarningInFunction
                    << "face[" << facei << "] = " << f
                    << " does not have three unique vertices" << endl;
            }
        }
    }

    // Detect doubled faces
    // do not touch the faces
    const labelListList& fFaces = this->faceFaces();
    newFacei = 0;
    forAll(faceLst, facei)
    {
        // skip already collapsed faces:
        if (faceMap[facei] < 0)
        {
            continue;
        }

        const Face& f = faceLst[facei];

        // duplicate face check
        bool okay = true;
        const labelList& neighbours = fFaces[facei];

        // Check if faceNeighbours use same points as this face.
        // Note: discards normal information - sides of baffle are merged.
        forAll(neighbours, neighI)
        {
            const label neiFacei = neighbours[neighI];

            if (neiFacei <= facei || faceMap[neiFacei] < 0)
            {
                // lower numbered faces already checked
                // skip neighbours that are themselves collapsed
                continue;
            }

            const Face& nei = faceLst[neiFacei];

            if (f == nei)
            {
                okay = false;

                if (verbose)
                {
                    WarningInFunction
                        << "faces share the same vertices:" << nl
                        << "    face[" << facei << "] : " << f << nl
                        << "    face[" << neiFacei << "] : " << nei << endl;
                    // printFace(Warning, "    ", f, points());
                    // printFace(Warning, "    ", nei, points());
                }

                break;
            }
        }

        if (okay)
        {
            faceMap[facei] = facei;
            ++newFacei;
        }
        else
        {
            faceMap[facei] = -1;
        }
    }

    // Phase 1: pack
    // Done to keep numbering constant in phase 1

    if (changed || newFacei < faceLst.size())
    {
        changed = true;

        if (verbose)
        {
            WarningInFunction
                << "Removed " << faceLst.size() - newFacei
                << " illegal faces." << endl;
        }

        // compress the face list
        newFacei = 0;
        forAll(faceLst, facei)
        {
            if (faceMap[facei] >= 0)
            {
                if (newFacei != facei)
                {
                    faceLst[newFacei] = faceLst[facei];
                }
                faceMap[newFacei] = facei;
                ++newFacei;
            }
        }

        faceLst.setSize(newFacei);
        remapFaces(faceMap);
    }
    faceMap.clear();

    // Topology can change because of renumbering
    ParentType::clearOut();
    return changed;
}


template<class Face>
Foam::label Foam::MeshedSurface<Face>::nTriangles() const
{
    if (faceTraits<Face>::isTri())
    {
        return ParentType::size();
    }

    return nTriangles
    (
        const_cast<List<label>&>(List<label>::null())
    );
}


template<class Face>
Foam::label Foam::MeshedSurface<Face>::nTriangles
(
    List<label>& faceMap
) const
{
    label nTri = 0;
    const List<Face>& faceLst = surfFaces();

    // Count triangles needed
    forAll(faceLst, facei)
    {
        nTri += faceLst[facei].nTriangles();
    }

    // Nothing to do
    if (nTri <= faceLst.size())
    {
        if (notNull(faceMap))
        {
            faceMap.clear();
        }
    }
    else if (notNull(faceMap))
    {
        // face map requested
        faceMap.setSize(nTri);

        nTri = 0;
        forAll(faceLst, facei)
        {
            label n = faceLst[facei].nTriangles();
            while (n-- > 0)
            {
                faceMap[nTri++] = facei;
            }
        }

        faceMap.setSize(nTri);
    }

    return nTri;
}


template<class Face>
Foam::label Foam::MeshedSurface<Face>::triangulate()
{
    if (faceTraits<Face>::isTri())
    {
        // Inplace triangulation of triFace/labelledTri surface = no-op
        return 0;
    }
    else
    {
        return triangulate
        (
            const_cast<List<label>&>(List<label>::null())
        );
    }
}


template<class Face>
Foam::label Foam::MeshedSurface<Face>::triangulate
(
    List<label>& faceMapOut
)
{
    if (faceTraits<Face>::isTri())
    {
        // Inplace triangulation of triFace/labelledTri surface = no-op
        if (notNull(faceMapOut))
        {
            faceMapOut.clear();
        }

        return 0;
    }

    label nTri = 0;
    label maxTri = 0;  // the maximum number of triangles for any single face
    List<Face>& faceLst = this->storedFaces();

    // determine how many triangles will be needed
    forAll(faceLst, facei)
    {
        const label n = faceLst[facei].nTriangles();
        if (maxTri < n)
        {
            maxTri = n;
        }
        nTri += n;
    }

    // nothing to do
    if (nTri <= faceLst.size())
    {
        if (notNull(faceMapOut))
        {
            faceMapOut.clear();
        }
        return 0;
    }

    List<Face>  newFaces(nTri);
    List<label> faceMap;

    // reuse storage from optional faceMap
    if (notNull(faceMapOut))
    {
        faceMap.transfer(faceMapOut);
    }
    faceMap.setSize(nTri);

    if (this->points().empty())
    {
        // triangulate without points
        // simple face triangulation around f[0]
        nTri = 0;
        forAll(faceLst, facei)
        {
            const Face& f = faceLst[facei];

            for (label fp = 1; fp < f.size() - 1; ++fp)
            {
                const label fp1 = f.fcIndex(fp);

                newFaces[nTri] = Face{f[0], f[fp], f[fp1]};
                faceMap[nTri] = facei;
                ++nTri;
            }
        }
    }
    else
    {
        // triangulate with points
        List<face> tmpTri(maxTri);

        nTri = 0;
        forAll(faceLst, facei)
        {
            // 'face' not '<Face>'
            const face& f = faceLst[facei];

            label nTmp = 0;
            f.triangles(this->points(), nTmp, tmpTri);
            for (label triI = 0; triI < nTmp; triI++)
            {
                newFaces[nTri] = Face
                (
                    static_cast<labelUList&>(tmpTri[triI])
                );
                faceMap[nTri] = facei;
                ++nTri;
            }
        }
    }

    // The number of *additional* faces
    nTri -= faceLst.size();

    faceLst.transfer(newFaces);
    remapFaces(faceMap);

    // optionally return the faceMap
    if (notNull(faceMapOut))
    {
        faceMapOut.transfer(faceMap);
    }
    faceMap.clear();

    // Topology can change because of renumbering
    ParentType::clearOut();

    return nTri;
}


template<class Face>
Foam::MeshedSurface<Face> Foam::MeshedSurface<Face>::subsetMesh
(
    const labelHashSet& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const pointField& locPoints = this->localPoints();
    const List<Face>& locFaces  = this->localFaces();


    // Fill pointMap, faceMap
    PatchTools::subsetMap(*this, include, pointMap, faceMap);

    // Create compact coordinate list and forward mapping array
    pointField newPoints(pointMap.size());
    labelList oldToNew(locPoints.size());
    forAll(pointMap, pointi)
    {
        newPoints[pointi] = locPoints[pointMap[pointi]];
        oldToNew[pointMap[pointi]] = pointi;
    }

    // create/copy a new zones list, each zone with zero size
    surfZoneList newZones(this->surfZones());
    forAll(newZones, zoneI)
    {
        newZones[zoneI].size() = 0;
    }

    // Renumber face node labels
    List<Face> newFaces(faceMap.size());
    forAll(faceMap, facei)
    {
        const label origFacei = faceMap[facei];
        newFaces[facei] = Face(locFaces[origFacei]);

        // Renumber labels for face
        Face& f = newFaces[facei];
        forAll(f, fp)
        {
            f[fp] = oldToNew[f[fp]];
        }
    }
    oldToNew.clear();

    // recalculate the zones start/size
    label newFacei = 0;
    label origEndI = 0;

    // adjust zone sizes
    forAll(newZones, zoneI)
    {
        surfZone& zone = newZones[zoneI];

        // adjust zone start
        zone.start() = newFacei;
        origEndI += zone.size();

        for (label facei = newFacei; facei < faceMap.size(); ++facei)
        {
            if (faceMap[facei] < origEndI)
            {
                ++newFacei;
            }
            else
            {
                break;
            }
        }

        // adjust zone size
        zone.size() = newFacei - zone.start();
    }


    // construct a sub-surface
    return MeshedSurface
    (
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newZones)
    );
}


template<class Face>
Foam::MeshedSurface<Face> Foam::MeshedSurface<Face>::subsetMesh
(
    const labelHashSet& include
) const
{
    labelList pointMap, faceMap;
    return subsetMesh(include, pointMap, faceMap);
}



template<class Face>
void Foam::MeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    ParentType::clearOut();

    this->storedPoints().transfer(surf.storedPoints());
    this->storedFaces().transfer(surf.storedFaces());
    this->storedZones().transfer(surf.storedZones());

    surf.clear();
}


template<class Face>
void Foam::MeshedSurface<Face>::transfer
(
    UnsortedMeshedSurface<Face>& surf
)
{
    clear();

    labelList faceMap;
    surfZoneList zoneLst = surf.sortedZones(faceMap);

    if (zoneLst.size() <= 1)
    {
        reset
        (
            xferMove(surf.storedPoints()),
            xferMove(surf.storedFaces()),
            Xfer<surfZoneList>()
        );
    }
    else
    {
        List<Face>& oldFaces = surf.storedFaces();
        List<Face> newFaces(faceMap.size());

        forAll(faceMap, facei)
        {
            newFaces[faceMap[facei]].transfer(oldFaces[facei]);
        }

        reset
        (
            xferMove(surf.storedPoints()),
            xferMove(newFaces),
            xferMove(zoneLst)
        );
    }

    faceMap.clear();
    surf.clear();
}


template<class Face>
Foam::Xfer<Foam::MeshedSurface<Face>>
Foam::MeshedSurface<Face>::xfer()
{
    return xferMove(*this);
}


template<class Face>
Foam::Xfer<Foam::List<Face>>
Foam::MeshedSurface<Face>::xferFaces()
{
    // Topology changed because of transfer
    ParentType::clearOut();

    return this->storedFaces().xfer();
}


template<class Face>
Foam::Xfer<Foam::List<Foam::point>>
Foam::MeshedSurface<Face>::xferPoints()
{
    // Topology changed because of transfer
    ParentType::clearOut();

    return this->storedPoints().xfer();
}


template<class Face>
Foam::Xfer<Foam::surfZoneList>
Foam::MeshedSurface<Face>::xferZones()
{
    return this->storedZones().xfer();
}


template<class Face>
bool Foam::MeshedSurface<Face>::read(const fileName& name)
{
    const word ext(name.ext());
    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();
        return read(unzipName, unzipName.ext());
    }

    return read(name, ext);
}


// Read from file in given format
template<class Face>
bool Foam::MeshedSurface<Face>::read
(
    const fileName& name,
    const word& ext
)
{
    clear();

    // read via selector mechanism
    transfer(New(name, ext)());
    return true;
}


template<class Face>
void Foam::MeshedSurface<Face>::write
(
    const Time& t,
    const word& surfName
) const
{
    MeshedSurfaceProxy<Face>(*this).write(t, surfName);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::operator=(const MeshedSurface& surf)
{
    clear();

    this->storedPoints() = surf.points();
    this->storedFaces()  = surf.surfFaces();
    this->storedZones()  = surf.surfZones();
}


template<class Face>
Foam::MeshedSurface<Face>::operator Foam::MeshedSurfaceProxy<Face>() const
{
    return MeshedSurfaceProxy<Face>
    (
        this->points(),
        this->surfFaces(),
        this->surfZones()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MeshedSurfaceZones.C"
#include "MeshedSurfaceIO.C"
#include "MeshedSurfaceNew.C"

// ************************************************************************* //

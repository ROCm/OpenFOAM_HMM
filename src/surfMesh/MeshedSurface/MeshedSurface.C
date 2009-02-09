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

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
#include "ListOps.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
// #include "surfMesh.H"
#include "primitivePatch.H"
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
    // handle 'native' format directly
    if (isNative(ext))
    {
        return true;
    }
    else
    {
        return checkSupport
        (
            readTypes() | SiblingType::readTypes(),
            ext,
            verbose,
            "reading"
        );
    }
}


template<class Face>
bool Foam::MeshedSurface<Face>::canWriteType
(
    const word& ext,
    const bool verbose
)
{
    // handle 'native' format directly
    if (isNative(ext))
    {
        return true;
    }

    return checkSupport(writeTypes(), ext, verbose, "writing");
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
    const MeshedSurface& surf
)
{
    if (debug)
    {
        Info<< "MeshedSurface::write"
            "(const fileName&, const MeshedSurface&) : "
            "writing to " << name
            << endl;
    }

    word ext = name.ext();

    // handle 'native' format directly
    if (isNative(ext))
    {
        surf.write(OFstream(name)());
        return;
    }

    typename writefileExtensionMemberFunctionTable::iterator mfIter =
        writefileExtensionMemberFunctionTablePtr_->find(ext);

    if (mfIter == writefileExtensionMemberFunctionTablePtr_->end())
    {
        FatalErrorIn
        (
            "MeshedSurface::write(const fileName&)"
        )   << "Unknown file extension " << ext << nl << nl
            << "Valid types are :" << endl
            << writeTypes()
            << exit(FatalError);
    }

    mfIter()(name, surf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface()
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst,
    const Xfer<surfRegionList>& regionLst
)
:
    ParentType(pointLst, faceLst),
    regions_(regionLst)
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst,
    const UList<label>& regionSizes,
    const UList<word>& regionNames
)
:
    ParentType(pointLst, faceLst)
{
    if (&regionSizes)
    {
        if (&regionNames)
        {
            addRegions(regionSizes, regionNames);
        }
        else
        {
            addRegions(regionSizes);
        }
    }
    else
    {
        oneRegion();
    }
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const polyBoundaryMesh& bMesh,
    const bool useGlobalPoints
)
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

    // use global/local points
    const pointField& bPoints =
    (
        useGlobalPoints ? mesh.points() : allBoundary.localPoints()
    );

    // global or local face addressing
    const List<Face>& bFaces =
    (
        useGlobalPoints ? allBoundary : allBoundary.localFaces()
    );


    // create region list
    surfRegionList newRegions(bPatches.size());

    label startFaceI = 0;
    forAll(bPatches, patchI)
    {
        const polyPatch& p = bPatches[patchI];

        newRegions[patchI] = surfRegion
        (
            p.name(),
            p.size(),
            startFaceI,
            patchI
        );

        startFaceI += p.size();
    }

    // must create with the same face type as the polyBoundaryMesh
    MeshedSurface<face> surf
    (
        xferCopy(bPoints),
        xferCopy(bFaces),
        xferMove(newRegions)
    );


    // must triangulate?
    if (this->isTri())
    {
        surf.triangulate();
        this->storedPoints().transfer(surf.storedPoints());

        // transcribe from face -> triFace
        List<face>&    origFaces = surf.storedFaces();
        List<triFace>  newFaces(origFaces.size());
        forAll(origFaces, faceI)
        {
            newFaces[faceI] = Face
            (
                static_cast<const UList<label>&>(origFaces[faceI])
            );
        }
        origFaces.clear();

        this->storedFaces().transfer(newFaces);
        regions_.transfer(surf.regions_);
    }
    else
    {
        this->transfer(surf);
    }
}


#if 0
// in preparation
Foam::MeshedSurface<Face>::MeshedSurface
(
    const surfMesh& sMesh
)
:
    ParentType(xferCopy(sMesh.points()), xferCopy(sMesh.faces()))
{
    const surfPatchList& sPatches = sMesh.boundaryMesh();

    // create regions list
    List<surfRegion> newRegions(sPatches.size());

    label startFaceI = 0;
    forAll(sPatches, patchI)
    {
        const surfPatch& p = sPatches[patchI];

        newRegions[patchI] = surfRegion
        (
            p.name(),
            p.size(),
            startFaceI,
            patchI
        );

        startFaceI += p.size();
    }

    regions_.transfer(newRegions);
}
#endif


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const UnsortedMeshedSurface<Face>& surf
)
{
    labelList faceMap;
    surfRegionList regionLst = surf.sortedRegions(faceMap);
    regions_.transfer(regionLst);

    const List<Face>& origFaces = surf.faces();
    List<Face> newFaces(origFaces.size());

    // this is somewhat like ListOps reorder and/or IndirectList
    forAll(newFaces, faceI)
    {
        newFaces[faceI] = origFaces[faceMap[faceI]];
    }

    reset(xferCopy(surf.points()), xferMove(newFaces));
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const fileName& name,
    const word& ext
)
{
    read(name, ext);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const fileName& name)
{
    read(name);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(Istream& is)
{
    read(is);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const Time& d)
{
    read(IFstream(findMeshName(d))());
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const MeshedSurface& surf)
:
    ParentType(surf),
    regions_(surf.regions_)
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<UnsortedMeshedSurface<Face> >& surf
)
{
    transfer(surf());
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const Xfer<MeshedSurface>& surf)
{
    transfer(surf());
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurface<Face>::~MeshedSurface()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::oneRegion(const word& name)
{
    word regionName(name);
    if (regionName.empty())
    {
        if (regions_.size())
        {
            regionName = regions_[0].name();
        }
        if (regionName.empty())
        {
            regionName = "region0";
        }
    }

    // set single default region
    regions_.setSize(1);
    regions_[0] = surfRegion
    (
        regionName,
        size(),         // region size
        0,              // region start
        0               // region index
    );
}


template<class Face>
void Foam::MeshedSurface<Face>::checkRegions()
{
    // extra safety, ensure we have at some regions
    // and they cover all the faces - fix start silently
    if (regions_.size() <= 1)
    {
        oneRegion();
    }
    else
    {
        label count = 0;
        forAll(regions_, regionI)
        {
            regions_[regionI].start() = count;
            count += regions_[regionI].size();
        }

        if (count < size())
        {
            WarningIn
            (
                "MeshedSurface::checkRegions()\n"
            )
                << "more face " << size() << " than regions " << count
                << " ... extending final region"
                << endl;

            regions_[regions_.size()-1].size() += count - size();
        }
        else if (count > size())
        {
            FatalErrorIn
            (
                "MeshedSurface::checkRegions()\n"
            )
                << "more regions " << count << " than faces " << size()
                << exit(FatalError);
        }
    }
}


template<class Face>
void Foam::MeshedSurface<Face>::sortFacesAndStore
(
    const Xfer<List<Face> >& unsortedFaces,
    const Xfer<List<label> >& regionIds,
    const bool sorted
)
{
    List<Face>  oldFaces(unsortedFaces);
    List<label> regions(regionIds);

    if (sorted)
    {
        // already sorted - simply transfer faces
        this->storedFaces().transfer(oldFaces);
    }
    else
    {
        // unsorted - determine the sorted order:
        // avoid SortableList since we discard the main list anyhow
        List<label> faceMap;
        sortedOrder(regions, faceMap);
        regions.clear();

        // sorted faces
        List<Face> newFaces(faceMap.size());
        forAll(faceMap, faceI)
        {
            // use transfer to recover memory if possible
            newFaces[faceI].transfer(oldFaces[faceMap[faceI]]);
        }
        this->storedFaces().transfer(newFaces);

    }
    regions.clear();
}


template<class Face>
void Foam::MeshedSurface<Face>::remapFaces
(
    const UList<label>& faceMap
)
{
    // recalculate the region start/size
    if (&faceMap && faceMap.size())
    {
        if (regions_.empty())
        {
            oneRegion();
        }
        else if (regions_.size() == 1)
        {
            // optimized for single region case
            regions_[0].size() = faceMap.size();
        }
        else
        {
            label newFaceI = 0;
            label origEndI = 0;
            forAll(regions_, regionI)
            {
                surfRegion& reg = regions_[regionI];

                // adjust region start
                reg.start() = newFaceI;
                origEndI += reg.size();

                for (label faceI = newFaceI; faceI < faceMap.size(); ++faceI)
                {
                    if (faceMap[faceI] < origEndI)
                    {
                        ++newFaceI;
                    }
                    else
                    {
                        break;
                    }
                }

                // adjust region size
                reg.size() = newFaceI - reg.start();
            }
        }
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::clear()
{
    ParentType::clear();
    regions_.clear();
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
    forAll(pointMap, pointI)
    {
        newPoints[pointI] = locPoints[pointMap[pointI]];
        oldToNew[pointMap[pointI]] = pointI;
    }

    // create/copy a new region list, each region with zero size
    surfRegionList newRegions(regions_);
    forAll(newRegions, regionI)
    {
        newRegions[regionI].size() = 0;
    }

    // Renumber face node labels
    List<Face> newFaces(faceMap.size());
    forAll(faceMap, faceI)
    {
        const label origFaceI = faceMap[faceI];
        newFaces[faceI] = Face(locFaces[origFaceI]);

        // Renumber labels for face
        Face& f = newFaces[faceI];
        forAll(f, fp)
        {
            f[fp] = oldToNew[f[fp]];
        }
    }
    oldToNew.clear();

    // recalculate the region start/size
    label newFaceI = 0;
    label origEndI = 0;

    // adjust region sizes
    forAll(newRegions, regionI)
    {
        surfRegion& reg = newRegions[regionI];

        // adjust region start
        reg.start() = newFaceI;
        origEndI += reg.size();

        for (label faceI = newFaceI; faceI < faceMap.size(); ++faceI)
        {
            if (faceMap[faceI] < origEndI)
            {
                ++newFaceI;
            }
            else
            {
                break;
            }
        }

        // adjust region size
        reg.size() = newFaceI - reg.start();
    }


    // construct a sub-surface
    return MeshedSurface
    (
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newRegions)
    );
}


template<class Face>
Foam::MeshedSurface<Face>
Foam::MeshedSurface<Face>::subsetMesh
(
    const labelHashSet& include
) const
{
    labelList pointMap, faceMap;
    return subsetMesh(include, pointMap, faceMap);
}


template<class Face>
void Foam::MeshedSurface<Face>::addRegions
(
    const UList<surfRegion>& regions,
    const bool cullEmpty
)
{
    label nRegion = 0;

    regions_.setSize(regions.size());
    forAll(regions_, regionI)
    {
        if (regions[regionI].size() || !cullEmpty)
        {
            regions_[nRegion] = surfRegion(regions[regionI], nRegion);
            nRegion++;
        }
    }
    regions_.setSize(nRegion);
}


template<class Face>
void Foam::MeshedSurface<Face>::addRegions
(
    const UList<label>& sizes,
    const UList<word>& names,
    const bool cullEmpty
)
{
    label start   = 0;
    label nRegion = 0;

    regions_.setSize(sizes.size());
    forAll(regions_, regionI)
    {
        if (sizes[regionI] || !cullEmpty)
        {
            regions_[nRegion] = surfRegion
            (
                names[regionI],
                sizes[regionI],
                start,
                nRegion
            );
            start += sizes[regionI];
            nRegion++;
        }
    }
    regions_.setSize(nRegion);
}


template<class Face>
void Foam::MeshedSurface<Face>::addRegions
(
    const UList<label>& sizes,
    const bool cullEmpty
)
{
    label start   = 0;
    label nRegion = 0;

    regions_.setSize(sizes.size());
    forAll(regions_, regionI)
    {
        if (sizes[regionI] || !cullEmpty)
        {
            regions_[nRegion] = surfRegion
            (
                word("region") + ::Foam::name(nRegion),
                sizes[regionI],
                start,
                nRegion
            );
            start += sizes[regionI];
            nRegion++;
        }
    }
    regions_.setSize(nRegion);
}


template<class Face>
void Foam::MeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    reset(xferMove(surf.storedPoints()), xferMove(surf.storedFaces()));
    regions_.transfer(surf.regions_);

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
    surfRegionList regionLst = surf.sortedRegions(faceMap);
    List<Face>& oldFaces = surf.storedFaces();

    List<Face> newFaces(faceMap.size());
    forAll(faceMap, faceI)
    {
        newFaces[faceI].transfer(oldFaces[faceMap[faceI]]);
    }
    faceMap.clear();

    reset(xferMove(surf.storedPoints()), xferMove(newFaces));
    regions_.transfer(regionLst);

    surf.clear();
}


// Read from file, determine format from extension
template<class Face>
bool Foam::MeshedSurface<Face>::read(const fileName& name)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();
        return read(unzipName, unzipName.ext());
    }
    else
    {
        return read(name, ext);
    }
}


// Read from file in given format
template<class Face>
bool Foam::MeshedSurface<Face>::read
(
    const fileName& name,
    const word& ext
)
{
    // handle 'native' format directly
    if (isNative(ext))
    {
        return read(IFstream(name)());
    }
    else
    {
        // use selector mechanism
        transfer(New(name, ext)());
        return true;
    }
}


template<class Face>
void Foam::MeshedSurface<Face>::write(const Time& d) const
{
    write(OFstream(findMeshName(d))());
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::operator=(const MeshedSurface& surf)
{
    clear();

    this->storedPoints() = surf.points();
    this->storedFaces()  = surf.faces();
    regions_ = surf.regions_;
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MeshedSurfaceIO.C"
#include "MeshedSurfaceNew.C"

// ************************************************************************* //

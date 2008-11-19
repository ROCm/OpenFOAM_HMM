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

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
#include "SortableList.H"
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
    wordHashSet known(2*fileExtensionConstructorTablePtr_->size());

    forAllIter
    (
        typename fileExtensionConstructorTable::iterator,
        *fileExtensionConstructorTablePtr_,
        iter
    )
    {
        known.insert(iter.key());
    }

    return known;
}


template<class Face>
Foam::wordHashSet Foam::MeshedSurface<Face>::writeTypes()
{
    wordHashSet supported(2*writefileExtensionMemberFunctionTablePtr_->size());

    forAllIter
    (
        typename writefileExtensionMemberFunctionTable::iterator,
        *writefileExtensionMemberFunctionTablePtr_,
        iter
    )
    {
        supported.insert(iter.key());
    }

    return supported;
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

    wordHashSet available = readTypes();
    available += SiblingType::readTypes();

    return checkSupport(available, ext, verbose, "reading");
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
    const fileName& fName,
    const bool verbose
)
{
    word ext = fName.ext();
    if (ext == "gz")
    {
        ext = fName.lessExt().ext();
    }
    return canReadType(ext, verbose);
}


template<class Face>
void Foam::MeshedSurface<Face>::write
(
    const fileName& fName,
    const MeshedSurface& surf
)
{
    if (debug)
    {
        Info<< "MeshedSurface::write(const fileName&, const MeshedSurface&) : "
               "writing MeshedSurface to " << fName
            << endl;
    }

    word ext = fName.ext();

    // handle 'native' format directly
    if (isNative(ext))
    {
        surf.write(OFstream(fName)());
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

    mfIter()(fName, surf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface()
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst,
    const xfer<surfGroupList>& patchLst
)
:
    ParentType(pointLst, faceLst),
    patches_(patchLst)
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst,
    const UList<label>& patchSizes,
    const UList<word>& patchNames,
    const UList<word>& patchTypes
)
:
    ParentType(pointLst, faceLst)
{
    surfGroupList newPatches(patchSizes.size());

    label start = 0;
    forAll(newPatches, patchI)
    {
        newPatches[patchI] = surfGroup
        (
            patchNames[patchI],
            patchSizes[patchI],
            start,
            patchI
        );

        start += patchSizes[patchI];
    }

    patches_.transfer(newPatches);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst,
    const UList<label>& regionIds,
    const Map<word>& regionNames
)
:
    ParentType(pointLst, faceLst)
{
    if (&regionIds && regionIds.size() && regionIds.size() != size())
    {
        FatalErrorIn
        (
            "MeshedSurface::MeshedSurface(\n"
            "(\n"
            "    const xfer<pointField>&,\n"
            "    const xfer<List<Face> >&,\n"
            "    const UList<label>& regionIds,\n"
            "    const Map<word>& regionNames\n"
            " )\n"
        )
            << "size mismatch : region and face sizes"
            << exit(FatalError);
    }
    else
    {
        sortFacesByRegion(regionIds, regionNames);
    }
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst,
    const UList<label>& regionIds,
    const HashTable<label>& nameToRegionMapping
)
:
    ParentType(pointLst, faceLst)
{
    if (&regionIds && regionIds.size() && regionIds.size() != size())
    {
        FatalErrorIn
        (
            "MeshedSurface::MeshedSurface(\n"
            "(\n"
            "    const xfer<pointField>&,\n"
            "    const xfer<List<Face> >&,\n"
            "    const UList<label>& regionIds,\n"
            "    const HashTable<label>& nameToRegionMapping\n"
            " )\n"
        )
            << "size mismatch : region and face sizes"
            << exit(FatalError);
    }
    else
    {
        Map<word> regionNames;
        forAllConstIter(HashTable<label>, nameToRegionMapping, iter)
        {
            regionNames.insert(iter(), iter.key());
        }

        sortFacesByRegion(regionIds, regionNames);
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


    // create patch list
    surfGroupList newPatches(bPatches.size());

    label startFaceI = 0;
    forAll(bPatches, patchI)
    {
        const polyPatch& p = bPatches[patchI];

        newPatches[patchI] = surfGroup
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
        xferMove(newPatches)
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
        patches_.transfer(surf.patches_);
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

    // create patch list
    List<surfGroup> newPatches(sPatches.size());

    label startFaceI = 0;
    forAll(sPatches, patchI)
    {
        const surfPatch& p = sPatches[patchI];

        newPatches[patchI] = surfGroup
        (
            p.name(),
            p.size(),
            startFaceI,
            patchI
        );

        startFaceI += p.size();
    }

    patches_.transfer(newPatches);
}
#endif


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const UnsortedMeshedSurface<Face>& surf
)
{
    labelList faceMap;
    surfGroupList patchLst = surf.sortedRegions(faceMap);
    patches_.transfer(patchLst);

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
    const fileName& fName,
    const word& ext
)
{
    read(fName, ext);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const fileName& fName)
{
    read(fName);
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
    // setDefaultPatches();
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const MeshedSurface& surf)
:
    ParentType(surf),
    patches_(surf.patches_)
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const xfer<UnsortedMeshedSurface<Face> >& surf
)
{
    transfer(surf());
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const xfer<MeshedSurface>& surf)
{
    transfer(surf());
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurface<Face>::~MeshedSurface()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::onePatch()
{
    // set single default patch
    patches_.setSize(1);
    patches_[0] = surfGroup
    (
        "patch0",
        size(),         // patch size
        0,              // patch start
        0               // patch index
    );
}


template<class Face>
void Foam::MeshedSurface<Face>::checkPatches()
{
    // extra safety, ensure we have at some patches,
    // and they cover all the faces
    // fix start silently
    if (patches_.size() > 1)
    {
        label count = 0;
        forAll(patches_, patchI)
        {
            patches_[patchI].start() = count;
            count += patches_[patchI].size();
        }

        if (count < size())
        {
            WarningIn
            (
                "MeshedSurface::checkPatches()\n"
            )
                << "more nFaces " << size()
                << " than patches " << count
                << " ... extending final patch"
                << endl;

            patches_[patches_.size()-1].size() += count - size();
        }
        else if (count > size())
        {
            FatalErrorIn
            (
                "MeshedSurface::checkPatches()\n"
            )
                << "more patches " << count
                << " than nFaces " << size()
                << exit(FatalError);
        }
    }
    else if (patches_.size() == 1)
    {
        // like onePatch, but preserve the name
        patches_[0].size() = size();
        patches_[0].start() = 0;
        if (!patches_[0].name().size())
        {
            patches_[0].name() = "patch0";
        }
    }
    else
    {
        onePatch();
    }
}


template<class Face>
void Foam::MeshedSurface<Face>::sortFacesByRegion
(
    const UList<label>& regionIds,
    const Map<word>& regionNames
)
{
    const List<Face>& unsortedFaces = this->faces();

    if (!&regionNames || !&regionIds || regionIds.size() == 0)
    {
        onePatch();
    }
    else if (regionIds.size() == unsortedFaces.size())
    {
        labelList faceMap;
        surfGroupList newPatches = UnsortedMeshedSurface<Face>::sortedRegions
        (
            regionIds,
            regionNames,
            faceMap
        );
        patches_.transfer(newPatches);

        // this is somewhat like ListOps reorder and/or IndirectList
        List<Face> newFaces(unsortedFaces.size());
        forAll(newFaces, faceI)
        {
            newFaces[faceI] = unsortedFaces[faceMap[faceI]];
        }
        faceMap.clear();

        this->storedFaces().transfer(newFaces);
    }
}


template<class Face>
void Foam::MeshedSurface<Face>::remapRegions(List<label>& faceMap)
{
    // recalculate the patch start/size
    if (faceMap.size())
    {
        label newFaceI = 0;
        label oldPatchEnd = 0;
        forAll(patches_, patchI)
        {
            surfGroup& p = patches_[patchI];

            // adjust patch start
            p.start() = newFaceI;
            oldPatchEnd += p.size();

            for (label faceI = newFaceI; faceI < faceMap.size(); ++faceI)
            {
                if (faceMap[faceI] < oldPatchEnd)
                {
                    ++newFaceI;
                }
                else
                {
                    break;
                }
            }

            // adjust patch size
            p.size() = newFaceI - p.start();
        }
        faceMap.clear();
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::clear()
{
    ParentType::clear();
    patches_.clear();
}


template<class Face>
Foam::MeshedSurface<Face> Foam::MeshedSurface<Face>::subsetMesh
(
    const UList<bool>& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const pointField& locPoints = this->localPoints();
    const List<Face>& locFaces  = this->localFaces();

    // Fill pointMap, faceMap
    this->subsetMap(include, pointMap, faceMap);

    // Create compact coordinate list and forward mapping array
    pointField newPoints(pointMap.size());
    labelList oldToNew(locPoints.size());
    forAll(pointMap, pointI)
    {
        newPoints[pointI] = locPoints[pointMap[pointI]];
        oldToNew[pointMap[pointI]] = pointI;
    }

    // create/copy a new patch list, each patch with zero size
    surfGroupList newPatches(patches_);
    forAll(newPatches, patchI)
    {
        newPatches[patchI].size() = 0;
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

    // recalculate the patch start/size
    label newFaceI = 0;
    label oldPatchEnd = 0;

    // adjust patch sizes
    forAll(newPatches, patchI)
    {
        surfGroup& p = newPatches[patchI];

        // adjust patch start
        p.start() = newFaceI;
        oldPatchEnd += p.size();

        for (label faceI = newFaceI; faceI < faceMap.size(); ++faceI)
        {
            if (faceMap[faceI] < oldPatchEnd)
            {
                ++newFaceI;
            }
            else
            {
                break;
            }
        }

        // adjust patch size
        p.size() = newFaceI - p.start();
    }


    // construct a sub-surface
    return MeshedSurface
    (
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newPatches)
    );
}


template<class Face>
Foam::MeshedSurface<Face> Foam::MeshedSurface<Face>::subsetMesh
(
    const UList<bool>& include
) const
{
    labelList pointMap, faceMap;
    return subsetMesh(include, pointMap, faceMap);
}


template<class Face>
void Foam::MeshedSurface<Face>::addPatches
(
    const UList<surfGroup>& patchLst
)
{
    patches_ = patchLst;
}


template<class Face>
void Foam::MeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    reset(xferMove(surf.storedPoints()), xferMove(surf.storedFaces()));
    patches_.transfer(surf.patches_);

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
    surfGroupList patchLst = surf.sortedRegions(faceMap);

    const List<Face>& oldFaces = surf.faces();
    List<Face>  newFaces(oldFaces.size());

    // this is somewhat like ListOps reorder and/or IndirectList
    forAll(newFaces, faceI)
    {
        newFaces[faceI] = oldFaces[faceMap[faceI]];
    }
    faceMap.clear();

    reset(xferMove(surf.storedPoints()), xferMove(newFaces));
    patches_.transfer(patchLst);

    surf.regions_.clear();
    surf.patches_.clear();
    surf.clear();
}


// Read from file, determine format from extension
template<class Face>
bool Foam::MeshedSurface<Face>::read(const fileName& fName)
{
    word ext = fName.ext();
    if (ext == "gz")
    {
        fileName unzipName = fName.lessExt();
        return read(unzipName, unzipName.ext());
    }
    else
    {
        return read(fName, ext);
    }
}


// Read from file in given format
template<class Face>
bool Foam::MeshedSurface<Face>::read
(
    const fileName& fName,
    const word& ext
)
{
    // handle 'native' format directly
    if (isNative(ext))
    {
        return read(IFstream(fName)());
    }
    else
    {
        // use selector mechanism
        transfer(New(fName, ext)());
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
    patches_ = surf.patches_;
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MeshedSurfaceCleanup.C"
#include "MeshedSurfaceIO.C"
#include "MeshedSurfaceNew.C"

// ************************************************************************* //

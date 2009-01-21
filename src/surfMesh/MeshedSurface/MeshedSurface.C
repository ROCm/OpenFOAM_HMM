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
    const Xfer<surfGroupList>& patchLst
)
:
    ParentType(pointLst, faceLst),
    patches_(patchLst)
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst,
    const UList<label>& patchSizes,
    const UList<word>& patchNames
)
:
    ParentType(pointLst, faceLst)
{
    if (&patchSizes)
    {
        if (&patchNames)
        {
            addPatches(patchSizes, patchNames);
        }
        else
        {
            addPatches(patchSizes);
        }
    }
    else
    {
        onePatch();
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
    patches_(surf.patches_)
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
void Foam::MeshedSurface<Face>::onePatch(const word& name)
{
    word patchName(name);
    if (patchName.empty())
    {
        if (patches_.size())
        {
            patchName = patches_[0].name();
        }
        if (patchName.empty())
        {
            patchName = "patch0";
        }
    }

    // set single default patch
    patches_.setSize(1);
    patches_[0] = surfGroup
    (
        patchName,
        size(),         // patch size
        0,              // patch start
        0               // patch index
    );
}


template<class Face>
void Foam::MeshedSurface<Face>::checkPatches()
{
    // extra safety, ensure we have at some patches,
    // and they cover all the faces - fix start silently
    if (patches_.size() <= 1)
    {
        onePatch();
    }
    else
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
                << "more face " << size() << " than patches " << count
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
                << "more patches " << count << " than faces " << size()
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
    // recalculate the patch start/size
    if (&faceMap && faceMap.size())
    {
        if (patches_.empty())
        {
            onePatch();
        }
        else if (patches_.size() == 1)
        {
            // optimized for one-patch case
            patches_[0].size() = faceMap.size();
        }
        else
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
        }
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
    const UList<surfGroup>& patches,
    const bool cullEmpty
)
{
    label nPatch = 0;

    patches_.setSize(patches.size());
    forAll(patches_, patchI)
    {
        if (patches[patchI].size() || !cullEmpty)
        {
            patches_[nPatch] = surfGroup(patches[patchI], nPatch);
            nPatch++;
        }
    }
    patches_.setSize(nPatch);
}


template<class Face>
void Foam::MeshedSurface<Face>::addPatches
(
    const UList<label>& sizes,
    const UList<word>& names,
    const bool cullEmpty
)
{
    label start  = 0;
    label nPatch = 0;

    patches_.setSize(sizes.size());
    forAll(patches_, patchI)
    {
        if (sizes[patchI] || !cullEmpty)
        {
            patches_[nPatch] = surfGroup
            (
                names[patchI],
                sizes[patchI],
                start,
                nPatch
            );
            start += sizes[patchI];
            nPatch++;
        }
    }
    patches_.setSize(nPatch);
}


template<class Face>
void Foam::MeshedSurface<Face>::addPatches
(
    const UList<label>& sizes,
    const bool cullEmpty
)
{
    label start  = 0;
    label nPatch = 0;

    patches_.setSize(sizes.size());
    forAll(patches_, patchI)
    {
        if (sizes[patchI] || !cullEmpty)
        {
            patches_[nPatch] = surfGroup
            (
                word("patch") + ::Foam::name(nPatch),
                sizes[patchI],
                start,
                nPatch
            );
            start += sizes[patchI];
            nPatch++;
        }
    }
    patches_.setSize(nPatch);
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
    List<Face>& oldFaces = surf.storedFaces();

    List<Face> newFaces(faceMap.size());
    forAll(faceMap, faceI)
    {
        newFaces[faceI].transfer(oldFaces[faceMap[faceI]]);
    }
    faceMap.clear();

    reset(xferMove(surf.storedPoints()), xferMove(newFaces));
    patches_.transfer(patchLst);

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
    patches_ = surf.patches_;
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MeshedSurfaceIO.C"
#include "MeshedSurfaceNew.C"

// ************************************************************************* //

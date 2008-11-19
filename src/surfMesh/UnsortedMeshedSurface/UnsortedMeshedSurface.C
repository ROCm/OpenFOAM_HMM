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
#include "boundBox.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "primitivePatch.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
Foam::wordHashSet Foam::UnsortedMeshedSurface<Face>::readTypes()
{
    wordHashSet supported(2*fileExtensionConstructorTablePtr_->size());

    forAllIter
    (
        typename fileExtensionConstructorTable::iterator,
        *fileExtensionConstructorTablePtr_,
        iter
    )
    {
        supported.insert(iter.key());
    }

    return supported;
}


template<class Face>
Foam::wordHashSet Foam::UnsortedMeshedSurface<Face>::writeTypes()
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


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::canReadType
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
    available += SiblingType::readTypes();;

    return checkSupport(available, ext, verbose, "reading");
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::canWriteType
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
bool Foam::UnsortedMeshedSurface<Face>::canRead
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
void Foam::UnsortedMeshedSurface<Face>::write
(
    const fileName& fName,
    const UnsortedMeshedSurface<Face>& surf
)
{
    if (debug)
    {
        Info<< "UnsortedMeshedSurface::write(const fileName&, const UnsortedMeshedSurface&) : "
               "writing UnsortedMeshedSurface to " << fName
            << endl;
    }

    const word ext = fName.ext();

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
            "UnsortedMeshedSurface::write"
            "(const fileName&, const UnsortedMeshedSurface&)"
        )   << "Unknown file extension " << ext << nl << nl
            << "Valid types are :" << endl
            << writeTypes()
            << exit(FatalError);
    }

    mfIter()(fName, surf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface()
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst,
    const xfer<List<label> >& regionIds,
    const xfer<surfPatchIdentifierList>& patchLst
)
:
    ParentType(pointLst, faceLst),
    regions_(regionIds),
    patches_(patchLst)
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst,
    const xfer<List<label> >& regionIds,
    const Map<word>& regionNames
)
:
    ParentType(pointLst, faceLst),
    regions_(regionIds)
{
    if (&regionNames)
    {
        // set patch names from (id => name) mapping
        setPatches(regionNames);
    }
    else
    {
        // find highest region ID and set patch names automatically
        setPatches();
    }
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst,
    const xfer<List<label> >& regionIds,
    const HashTable<label>& labelToRegion
)
:
    ParentType(pointLst, faceLst),
    regions_(regionIds)
{
    // set patch names from (name => id) mapping
    setPatches(labelToRegion);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst
)
:
    ParentType(pointLst, faceLst)
{
    onePatch();
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const polyBoundaryMesh& bMesh,
    const bool useGlobalPoints
)
{
    // creating via MeshedSurface is the easiest
    MeshedSurface<Face> surf(bMesh, useGlobalPoints);
    transfer(surf);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const MeshedSurface<Face>& surf
)
:
    ParentType(xferCopy(surf.points()), xferCopy(surf.faces()))
{
    const surfGroupList& patchLst = surf.patches();

    regions_.setSize(size());
    patches_.setSize(patchLst.size());

    label faceI = 0;
    forAll(patchLst, patchI)
    {
        patches_[patchI] = patchLst[patchI];

        forAll(patchLst[patchI], patchFaceI)
        {
            regions_[faceI++] = patchI;
        }
    }
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const fileName& fName,
    const word& ext
)
{
    read(fName, ext);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface(const fileName& fName)
{
    read(fName);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface(Istream& is)
{
    read(is);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface(const Time& d)
{
    read(IFstream(findMeshName(d))());
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const UnsortedMeshedSurface<Face>& surf
)
:
    ParentType(xferCopy(surf.points()), xferCopy(surf.faces())),
    regions_(surf.regions_),
    patches_(surf.patches_)
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const xfer<UnsortedMeshedSurface<Face> >& surf
)
{
    transfer(surf());
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const xfer<MeshedSurface<Face> >& surf
)
{
    transfer(surf());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::UnsortedMeshedSurface<Face>::~UnsortedMeshedSurface()
{}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::onePatch()
{
    regions_.setSize(size());
    regions_ = 0;

    // set single default patch
    patches_.setSize(1);
    patches_[0] = PatchRegionType("patch0", 0);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setPatches(const label maxPatch)
{
    patches_.setSize(maxPatch+1);

    forAll(patches_, patchI)
    {
        patches_[patchI] = PatchRegionType
        (
            "patch" + ::Foam::name(patchI),
            patchI
        );
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setPatches()
{
    label maxPatch = 0;

    // find the max region that occurs
    forAll(regions_, faceI)
    {
        const label regId = regions_[faceI];

        if (maxPatch < regId)
        {
            maxPatch = regId;
        }
    }

    setPatches(maxPatch);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setPatches
(
    const Map<word>& regionNames,
    const label maxPatchHint
)
{
    label maxPatch = maxPatchHint;

    // determine max patch ID if required
    if (maxPatchHint < 0)
    {
        maxPatch = 0;
        forAllConstIter(Map<word>, regionNames, iter)
        {
            if (maxPatch < iter.key())
            {
                maxPatch = iter.key();
            }
        }
    }


    // Info<< "setPatches with maxPatch: " << maxPatch << endl;

    patches_.setSize(maxPatch+1);

    forAll(patches_, patchI)
    {
        Map<word>::const_iterator findPatch = regionNames.find(patchI);
        word patchName;

        if (findPatch != regionNames.end())
        {
            patchName = findPatch();
        }
        else
        {
            patchName = "patch" + ::Foam::name(patchI);
        }

        patches_[patchI] = PatchRegionType
        (
            patchName,
            patchI
        );
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setPatches
(
    const HashTable<label>& groupToPatch
)
{
    // determine max patch Id
    label maxPatch = 0;
    Map<word> regionNames;

    forAllConstIter(HashTable<label>, groupToPatch, iter)
    {
        regionNames.insert(iter(), iter.key());
        if (maxPatch < iter())
        {
            maxPatch = iter();
        }
    }

    setPatches(regionNames, maxPatch);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::remapRegions(List<label>& faceMap)
{
    // re-assign the region Ids
    if (faceMap.size())
    {
        forAll(faceMap, faceI)
        {
            faceMap[faceI] = regions_[faceMap[faceI]];
        }
        regions_.transfer(faceMap);
    }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setSize(const label s)
{
    ParentType::setSize(s);
    // if regions extend: set with last patchId
    regions_.setSize(s, patches_.size() - 1);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::clear()
{
    ParentType::clear();
    regions_.clear();
    patches_.clear();
}


template<class Face>
Foam::surfGroupList Foam::UnsortedMeshedSurface<Face>::sortedRegions
(
    labelList& faceMap
) const
{
    // supply some patch names
    Map<word> patchNames;
    forAll(patches_, patchI)
    {
        patchNames.insert(patchI, patches_[patchI].name());
    }

    return sortedPatchRegions(regions_, patchNames, faceMap);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face> Foam::UnsortedMeshedSurface<Face>::subsetMesh
(
    const UList<bool>& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const pointField&  locPoints = this->localPoints();
    const List<Face>&  locFaces  = this->localFaces();

    // Fill pointMap, faceMap
    this->subsetMap(include, pointMap, faceMap);

    // Create compact coordinate list and forward mapping array
    pointField newPoints(pointMap.size());
    labelList  oldToNew(locPoints.size());
    forAll(pointMap, pointI)
    {
        newPoints[pointI] = locPoints[pointMap[pointI]];
        oldToNew[pointMap[pointI]] = pointI;
    }

    // Renumber face node labels and compact
    List<Face>  newFaces(faceMap.size());
    List<label> newRegions(faceMap.size());

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

        newRegions[faceI] = regions_[origFaceI];
    }
    oldToNew.clear();

    // construct a sub-surface
    return UnsortedMeshedSurface
    (
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newRegions),
        xferCopy(patches_)
    );
}


template<class Face>
Foam::UnsortedMeshedSurface<Face> Foam::UnsortedMeshedSurface<Face>::subsetMesh
(
    const UList<bool>& include
) const
{
    labelList pointMap, faceMap;
    return subsetMesh(include, pointMap, faceMap);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::reset
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst,
    const xfer<List<label> >& regionIds
)
{
    ParentType::reset(pointLst, faceLst);

    if (&regionIds)
    {
        regions_.transfer(regionIds());
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::transfer
(
    UnsortedMeshedSurface<Face>& surf
)
{
    reset
    (
        xferMove(surf.storedPoints()),
        xferMove(surf.storedFaces()),
        xferMove(surf.regions_)
    );
    patches_.transfer(surf.patches_);

    surf.clear();
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    surfGroupList& patchLst = surf.patches_;

    reset(xferMove(surf.storedPoints()), xferMove(surf.storedFaces()));

    regions_.setSize(size());
    patches_.setSize(patchLst.size());

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        // copy info
        patches_[patchI] = patchLst[patchI];

        forAll(patchLst[patchI], patchFaceI)
        {
            regions_[faceIndex++] = patchI;
        }
    }

    patchLst.clear();
    surf.clear();
}


// Read from file, determine format from extension
template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::read(const fileName& fName)
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
bool Foam::UnsortedMeshedSurface<Face>::read
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
void Foam::UnsortedMeshedSurface<Face>::write(const Time& d) const
{
    write(OFstream(findMeshName(d))());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::operator=
(
    const UnsortedMeshedSurface<Face>& surf
)
{
    clear();

    this->storedPoints() = surf.points();
    this->storedFaces()  = surf.faces();
    regions_ = surf.regions_;
    patches_ = surf.patches_;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UnsortedMeshedSurfaceCleanup.C"
#include "UnsortedMeshedSurfaceIO.C"
#include "UnsortedMeshedSurfaceNew.C"

// ************************************************************************* //

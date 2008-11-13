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

// File extension for 'native' raw format
#undef  nativeSurfaceExt
#define nativeSurfaceExt "ofs"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
Foam::string Foam::UnsortedMeshedSurface<Face>::getLineNoComment(IFstream& is)
{
    string line;
    do
    {
        is.getLine(line);
    }
    while ((line.size() == 0 || line[0] == '#') && is.good());

    return line;
}

template<class Face>
Foam::fileName Foam::UnsortedMeshedSurface<Face>::triSurfInstance(const Time& d)
{
    fileName foamName(d.caseName() + "." + nativeSurfaceExt);

    // Search back through the time directories list to find the time
    // closest to and lower than current time

    instantList ts = d.times();
    label i;

    for (i=ts.size()-1; i>=0; i--)
    {
        if (ts[i].value() <= d.timeOutputValue())
        {
            break;
        }
    }

    // Noting that the current directory has already been searched
    // for mesh data, start searching from the previously stored time directory

    if (i>=0)
    {
        for (label j=i; j>=0; j--)
        {
            if (file(d.path()/ts[j].name()/typeName/foamName))
            {
                if (debug)
                {
                    Pout<< " UnsortedMeshedSurface::triSurfInstance(const Time& d)"
                        << "reading " << foamName
                        << " from " << ts[j].name()/typeName
                        << endl;
                }

                return ts[j].name();
            }
        }
    }

    if (debug)
    {
        Pout<< " UnsortedMeshedSurface::triSurfInstance(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }

    return "constant";
}


template<class Face>
Foam::fileName Foam::UnsortedMeshedSurface<Face>::triSurfName(const Time& d)
{
    fileName foamName(d.caseName() + "." + nativeSurfaceExt);

    // Search back through the time directories list to find the time
    // closest to and lower than current time

    instantList ts = d.times();
    label i;

    for (i=ts.size()-1; i>=0; i--)
    {
        if (ts[i].value() <= d.timeOutputValue())
        {
            break;
        }
    }

    // Noting that the current directory has already been searched
    // for mesh data, start searching from the previously stored time directory

    if (i>=0)
    {
        for (label j=i; j>=0; j--)
        {
            fileName testName(d.path()/ts[j].name()/typeName/foamName);

            if (file(testName))
            {
                if (debug)
                {
                    Pout<< " UnsortedMeshedSurface::triSurfName(const Time& d)"
                        << "reading " << foamName
                        << " from " << ts[j].name()/typeName
                        << endl;
                }

                return testName;
            }
        }
    }

    if (debug)
    {
        Pout<< " UnsortedMeshedSurface::triSurfName(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }

    return d.path()/"constant"/typeName/foamName;
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::canRead
(
    const word& ext,
    const bool verbose
)
{
    // perhaps sent an entire name
    word fExt(ext);

    string::size_type dot = fExt.find_last_of(".");
    if (dot != string::npos)
    {
        fExt = fExt.substr(dot+1);
    }

    // handle 'native' format directly
    if (fExt == nativeSurfaceExt)
    {
        return true;
    }

    typename fileExtensionConstructorTable::iterator cstrIter =
        fileExtensionConstructorTablePtr_->find(fExt);

    // would be nice to have information about which format this actually is
    if (cstrIter == fileExtensionConstructorTablePtr_->end())
    {
        if (verbose)
        {
            SortableList<word> known
            (
                fileExtensionConstructorTablePtr_->toc()
            );

            Info<<"Unknown file extension for reading: " << fExt << nl;
            // compact output:
            Info<<"Valid types: ( " << nativeSurfaceExt;
            forAll(known, i)
            {
                Info<<" " << known[i];
            }
            Info<<" )" << endl;
        }
        return false;
    }

    return true;
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::canWrite
(
    const word& ext,
    const bool verbose
)
{
    // perhaps sent an entire name
    word fExt(ext);

    string::size_type dot = ext.find_last_of(".");
    if (dot != string::npos)
    {
        fExt = ext.substr(dot+1);
    }

    // handle 'native' format directly
    if (fExt == nativeSurfaceExt)
    {
        return true;
    }

    typename writefileExtensionMemberFunctionTable::iterator mfIter =
        writefileExtensionMemberFunctionTablePtr_->find(fExt);

    if (mfIter == writefileExtensionMemberFunctionTablePtr_->end())
    {
        if (verbose)
        {
            SortableList<word> known
            (
                writefileExtensionMemberFunctionTablePtr_->toc()
            );

            Info<<"Unknown file extension for writing: " << fExt << nl;
            // compact output:
            Info<<"Valid types: ( " << nativeSurfaceExt;
            forAll(known, i)
            {
                Info<<" " << known[i];
            }
            Info<<" )" << endl;
        }

        return false;
    }

    return true;
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
    if (ext == nativeSurfaceExt)
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
            << writefileExtensionMemberFunctionTablePtr_->toc()
            << exit(FatalError);
    }

    mfIter()(fName, surf);
}


// Returns patch info.
// Sets faceMap to the indexing according to patch numbers.
// Patch numbers start at 0.
template<class Face>
Foam::surfGroupList Foam::UnsortedMeshedSurface<Face>::sortedRegions
(
    const UList<label>& regionLst,
    const Map<word>& patchNames,
    labelList& faceMap
)
{
    // determine sort order according to region numbers

    // std::sort() really seems to mix up the order.
    // and std::stable_sort() might take too long / too much memory

    // Assuming that we have relatively fewer regions compared to the
    // number of items, just do it ourselves

    // step 1: get region sizes and store (regionId => patchI)
    Map<label> regionLookup;
    forAll(regionLst, faceI)
    {
        const label regId = regionLst[faceI];

        Map<label>::iterator iter = regionLookup.find(regId);
        if (iter == regionLookup.end())
        {
            regionLookup.insert(regId, 1);
        }
        else
        {
            iter()++;
        }
    }

    // step 2: assign start/size (and name) to the newPatches
    // re-use the lookup to map (regionId => patchI)
    surfGroupList patchLst(regionLookup.size());
    label patchStart = 0;
    label patchI = 0;
    forAllIter(Map<label>, regionLookup, iter)
    {
        label regId = iter.key();

        word patchName;
        Map<word>::const_iterator iter2 = patchNames.find(regId);
        if (iter2 == patchNames.end())
        {
            patchName = word("patch") + ::Foam::name(patchI);
        }
        else
        {
            patchName = iter2();
        }

        patchLst[patchI] = surfGroup
        (
            patchName,
            0,           // initialize with zero size
            patchStart,
            patchI
        );

        // increment the start for the next patch
        // and save the (regionId => patchI) mapping
        patchStart += iter();
        iter() = patchI++;
    }


    // step 3: build the re-ordering
    faceMap.setSize(regionLst.size());

    forAll(regionLst, faceI)
    {
        label patchI = regionLookup[regionLst[faceI]];

        faceMap[faceI] = patchLst[patchI].start() + patchLst[patchI].size()++;
    }

    // with reordered faces registered in faceMap
    return patchLst;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface()
:
    ParentType(List<Face>(), pointField())
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
    ParentType(List<Face>(), pointField()),
    regions_(regionIds),
    patches_(patchLst)
{
    points().transfer(pointLst());
    faces().transfer(faceLst());
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst,
    const xfer<List<label> >& regionIds,
    const Map<word>& regionNames
)
:
    ParentType(List<Face>(), pointField()),
    regions_(regionIds)
{
    faces().transfer(faceLst());
    points().transfer(pointLst());

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
    ParentType(List<Face>(), pointField()),
    regions_(regionIds)
{
    points().transfer(pointLst());
    faces().transfer(faceLst());

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
    ParentType(List<Face>(), pointField()),
    regions_(faceLst().size(), 0),     // single default patch
    patches_()
{
    points().transfer(pointLst());
    faces().transfer(faceLst());

    setPatches(0);
}


#if 0
template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const polyBoundaryMesh& bMesh,
    const bool useGlobalPoints
)
:
    ParentType(List<Face>(), pointField())
{
    const polyMesh& mesh = bMesh.mesh();
    const polyPatchList& bPatches = bMesh;
    const label nIntFaces = mesh.nInternalFaces();

    List<PatchRegionType> newPatches(bPatches.size());

    // Get patch for all of outside
    primitivePatch allBoundary
    (
        SubList<Face>
        (
            mesh.faces(),
            mesh.nFaces() - nIntFaces,
            nIntFaces
        ),
        mesh.points()
    );

    List<Face>  newFaces(allBoundary.size());
    List<label> newRegions(allBoundary.size());

    if (useGlobalPoints)
    {
        // copy in the global points
        points() = mesh.points();
    }
    else
    {
        // copy in the local points
        points() = allBoundary.localPoints();
    }

    // global or local face addressing
    const List<Face>& bfaces =
    (
        useGlobalPoints
      ? allBoundary
      : allBoundary.localFaces()
    );

    label faceIndex = 0;
    forAll(bPatches, patchI)
    {
        const polyPatch& p = bPatches[patchI];

        newPatches[patchI] = PatchRegionType
        (
            bPatches[patchI].name(),
            patchI
        );

        forAll(p, patchFaceI)
        {
            newFaces[faceIndex] = bfaces[faceIndex];
            newRegions[faceIndex] = patchI;
            faceIndex++;
        }
    }

    faces().transfer(newFaces);
    regions_.transfer(newRegions);
    patches_.transfer(newPatches);
}
#endif


#if 0
template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const MeshedSurface<Face>& surf
)
:
    ParentType(List<Face>(), surf.points())
{
    const List<Face>&   origFaces = surf.faces();
    const surfGroupList& patchLst = surf.patches();

    List<Face>   newFaces(origFaces.size());
    List<label>  newRegions(origFaces.size());
    List<PatchRegionType> newPatches(patchLst.size());

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        newPatches[patchI] = patchLst[patchI];
        forAll(patchLst[patchI], patchFaceI)
        {
            newFaces[faceIndex] = origFaces[faceIndex];
            newRegions[faceIndex] = patchI;
            faceIndex++;
        }
    }

    faces().transfer(newFaces);
    regions_.transfer(newRegions);
    patches_.transfer(newPatches);
}
#endif

template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const fileName& fName,
    const word& ext
)
:
    ParentType(List<Face>(), pointField())
{
    read(fName, ext);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const fileName& fName
)
:
    ParentType(List<Face>(), pointField())
{
    read(fName, fName.ext());
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    Istream& is
)
:
    ParentType(List<Face>(), pointField())
{
    read(is);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface(const Time& d)
:
    ParentType(List<Face>(), pointField())
{
    read(IFstream(triSurfName(d))());
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const UnsortedMeshedSurface<Face>& surf
)
:
    ParentType(surf.faces(), surf.points()),
    regions_(surf.regions_),
    patches_(surf.patches_)
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const xfer<UnsortedMeshedSurface<Face> >& surf
)
:
    ParentType(List<Face>(), pointField()),
    regions_(),
    patches_()
{
    transfer(surf());
}

#if 0
template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const xfer<MeshedSurface<Face> >& surf
)
:
    ParentType(List<Face>(), pointField()),
    regions_(),
    patches_()
{
    transfer(surf());
}
#endif

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face>
Foam::UnsortedMeshedSurface<Face>::~UnsortedMeshedSurface()
{}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setSize(const label s)
{
    ParentType::setSize(s);
    regions_.setSize(s);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::clear()
{
    ParentType::clearOut();

    points().clear();
    faces().clear();
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

    return sortedRegions(regions_, patchNames, faceMap);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::movePoints(const pointField& newPoints)
{
    // Remove all geometry dependent data
    ParentType::clearTopology();

    // Adapt for new point position
    ParentType::movePoints(newPoints);

    // Copy new points
    points() = newPoints;
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::scalePoints(const scalar& scaleFactor)
{
    // avoid bad scaling
    if (scaleFactor > 0 && scaleFactor != 1.0)
    {
        // Remove all geometry dependent data
        ParentType::clearTopology();

        // Adapt for new point position
        ParentType::movePoints(pointField());

        points() *= scaleFactor;
    }
}


template<class Face>
Foam::UnsortedMeshedSurface<Face> Foam::UnsortedMeshedSurface<Face>::subsetMesh
(
    const UList<bool>& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const pointField& locPoints = ParentType::localPoints();
    const List<Face>& locFaces  = ParentType::localFaces();

    // Fill pointMap, faceMap
    ParentType::subsetMap(include, pointMap, faceMap);

    // Create compact coordinate list and forward mapping array
    pointField newPoints(pointMap.size());
    labelList oldToNew(locPoints.size());
    forAll(pointMap, pointI)
    {
        newPoints[pointI] = locPoints[pointMap[pointI]];
        oldToNew[pointMap[pointI]] = pointI;
    }

    // Renumber face node labels and compact
    List<Face> newFaces(faceMap.size());

    forAll(faceMap, faceI)
    {
        // Get old vertex labels
        const Face& oldFace = locFaces[faceMap[faceI]];

        newFaces[faceI] = Face(oldFace);

        // Renumber labels for face
        Face& f = newFaces[faceI];
        forAll(f, fp)
        {
            f[fp] = oldToNew[oldFace[fp]];
        }
    }

    // construct a sub-surface
    UnsortedMeshedSurface<Face> subSurf;
    (
        xferMove(newPoints),
        xferMove(newFaces),
        xferCopy(patches_)
    );

    return subSurf;
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::transfer
(
    UnsortedMeshedSurface<Face>& surf
)
{
    clear();

    faces().transfer(surf.faces());
    points().transfer(surf.points());
    regions_.transfer(surf.regions_);
    patches_.transfer(surf.patches_);

    surf.clear();
}


#if 0
template<class Face>
void Foam::UnsortedMeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    surfGroupList& patchLst = surf.patches();

    clear();

    points().transfer(surf.points());
    faces().transfer(surf.faces());
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
#endif


// Read from file in given format
template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::read
(
    const fileName& fName
)
{
    clear();
    const word ext = fName.ext();

    // handle 'native' format directly
    if (ext == nativeSurfaceExt)
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


// Read from file in given format
template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::read
(
    const fileName& fName,
    const word& ext
)
{
    clear();

    // handle 'native' format directly
    if (ext == nativeSurfaceExt)
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
    write(OFstream(triSurfName(d))());
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::writeStats(Ostream& os) const
{
//    os  << "Faces        : " << size() << endl
//        << "Edges        : " << nEdges() << endl
//        << "Vertices     : " << nPoints() << endl
//        << "Bounding Box : " << boundBox(localPoints(), false) << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::operator=
(
    const UnsortedMeshedSurface<Face>& surf
)
{
    clear();
    faces()  = surf.faces();
    points() = surf.points();
    regions_ = surf.regions_;
    patches_ = surf.patches_;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UnsortedMeshedSurfaceCleanup.C"
#include "UnsortedMeshedSurfaceIO.C"
#include "UnsortedMeshedSurfaceNew.C"

#undef nativeSurfaceExt

// ************************************************************************* //

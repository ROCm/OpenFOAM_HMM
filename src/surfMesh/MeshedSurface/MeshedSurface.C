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
#include "demandDrivenData.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
#include "boundBox.H"
#include "SortableList.H"
#include "ListOps.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
// #include "surfMesh.H"
#include "primitivePatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// File extension for 'native' raw format
#undef  nativeSurfaceExt
#define nativeSurfaceExt "ofs"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
Foam::fileName Foam::MeshedSurface<Face>::triSurfInstance(const Time& d)
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
                    Pout<< " MeshedSurface::triSurfInstance(const Time& d)"
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
        Pout<< " MeshedSurface::triSurfInstance(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }

    return "constant";
}


template<class Face>
Foam::fileName Foam::MeshedSurface<Face>::triSurfName(const Time& d)
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
                    Pout<< " MeshedSurface::triSurfName(const Time& d)"
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
        Pout<< " MeshedSurface::triSurfName(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }

    return d.path()/"constant"/typeName/foamName;
}

template<class Face>
bool Foam::MeshedSurface<Face>::canRead(const word& ext, const bool verbose)
{
    // handle 'native' format directly
    if (ext == nativeSurfaceExt)
    {
        return true;
    }

    return UnsortedMeshedSurface<Face>::canRead(ext, verbose);
}


template<class Face>
bool Foam::MeshedSurface<Face>::canWrite(const word& ext, const bool verbose)
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
            "MeshedSurface::write(const fileName&)"
        )   << "Unknown file extension " << ext << nl << nl
            << "Valid types are :" << endl
            << writefileExtensionMemberFunctionTablePtr_->toc()
            << exit(FatalError);
    }

    mfIter()(fName, surf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface()
:
    ParentType(List<Face>(), pointField())
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<Face> >& faceLst,
    const xfer<surfGroupList>& patchLst
)
:
    ParentType(List<Face>(), pointField()),
    patches_(patchLst)
{
    points().transfer(pointLst());
    faces().transfer(faceLst());
}


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
    ParentType(List<Face>(), pointField())
{
    points().transfer(pointLst());
    faces().transfer(faceLst());

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
    ParentType(List<Face>(), pointField())
{
    points().transfer(pointLst());
    faces().transfer(faceLst());

    if
    (
        &regionIds
     && regionIds.size() != 0
     && regionIds.size() != nFaces()
    )
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
            << "size mismatch : regionIds.size() != nFaces()"
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
    ParentType(List<Face>(), pointField())
{
    points().transfer(pointLst());
    faces().transfer(faceLst());

    if (regionIds.size() != nFaces())
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
            << "size mismatch : regionIds.size() != nFaces()"
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
:
    ParentType(List<Face>(), pointField())
{
    const polyMesh& mesh = bMesh.mesh();
    const polyPatchList& bPatches = bMesh;
    const label nIntFaces = mesh.nInternalFaces();

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

    if (useGlobalPoints)
    {
        // copy in the global points and the global face addressing
        points() = mesh.points();
        faces() = allBoundary;
    }
    else
    {
        // copy in the local points and the local face addressing
        points() = allBoundary.localPoints();
        faces() = allBoundary.localFaces();
    }

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

    patches_.transfer(newPatches);
}


#if 0
// in preparation
Foam::MeshedSurface<Face>::MeshedSurface
(
    const surfMesh& sMesh
)
:
    ParentType(List<Face>(sMesh.faces()), sMesh.points())
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
:
    ParentType(List<Face>(), surf.points())
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

    faces().transfer(newFaces);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
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
Foam::MeshedSurface<Face>::MeshedSurface
(
    const fileName& fName
)
:
    ParentType(List<Face>(), pointField())
{
    read(fName, fName.ext());
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(Istream& is)
:
    ParentType(List<Face>(), pointField())
{
    read(is);
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const Time& d)
:
    ParentType(List<Face>(), pointField())
{
    read(IFstream(triSurfName(d))());
    // setDefaultPatches();
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const MeshedSurface& surf)
:
    ParentType(surf.faces(), surf.points()),
    patches_(surf.patches_)
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const xfer<UnsortedMeshedSurface<Face> >& surf
)
:
    ParentType(List<Face>(), pointField())
{
    transfer(surf());
}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface(const xfer<MeshedSurface>& surf)
:
    ParentType(List<Face>(), pointField())
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
        patches_[0].start() = 0;
        patches_[0].size() = size();
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
    const List<Face>& unsortedFaces = faces();

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

        faces().transfer(newFaces);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::clear()
{
    ParentType::clearOut();

    points().clear();
    faces().clear();
    patches_.clear();
}


template<class Face>
void Foam::MeshedSurface<Face>::movePoints(const pointField& newPoints)
{
    // Remove all geometry dependent data
    ParentType::clearTopology();

    // Adapt for new point position
    ParentType::movePoints(newPoints);

    // Copy new points
    points() = newPoints;
}


template<class Face>
void Foam::MeshedSurface<Face>::scalePoints(const scalar& scaleFactor)
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
Foam::MeshedSurface<Face> Foam::MeshedSurface<Face>::subsetMesh
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

    // create a new patch list
    surfGroupList newPatches(patches_);
    forAll(newPatches, patchI)
    {
        newPatches[patchI].size() = 0;
    }

    // Renumber face node labels and compact
    List<Face> newFaces(faceMap.size());

    forAll(faceMap, faceI)
    {
        // Get old vertex labels
        label origFaceI = faceMap[faceI];
        const Face& oldFace = locFaces[origFaceI];

        newFaces[faceI] = Face(oldFace);

        // Renumber labels for face
        Face& f = newFaces[faceI];
        forAll(f, fp)
        {
            f[fp] = oldToNew[oldFace[fp]];
        }

        // adjust patch sizes
        forAllReverse (newPatches, patchI)
        {
            if
            (
                origFaceI >= patches_[patchI].start()
             && patches_[patchI].size()
            )
            {
                newPatches[patchI].size()++;
                break;
            }
        }
    }

    oldToNew.clear();

    // adjust patch start
    label startFaceI = 0;
    forAll(newPatches, patchI)
    {
        newPatches[patchI].start() = startFaceI;
        startFaceI += newPatches[patchI].size();
    }

    // construct a sub-surface
    MeshedSurface subSurf
    (
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newPatches)
    );

    return subSurf;
}


template<class Face>
void Foam::MeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    clear();

    points().transfer(surf.points());
    faces().transfer(surf.faces());
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
    points().transfer(surf.points());

    labelList faceMap;
    surfGroupList patchLst = surf.sortedRegions(faceMap);
    patches_.transfer(patchLst);

    surf.regions_.clear();
    surf.patches_.clear();

    List<Face>& oldFaces = surf.faces();
    List<Face>  newFaces(oldFaces.size());

    // this is somewhat like ListOps reorder and/or IndirectList
    forAll(newFaces, faceI)
    {
        newFaces[faceI] = oldFaces[faceMap[faceI]];
    }

    faces().transfer(newFaces);

    surf.clear();
}


// Read from file in given format
template<class Face>
bool Foam::MeshedSurface<Face>::read
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
bool Foam::MeshedSurface<Face>::read
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
void Foam::MeshedSurface<Face>::write(const Time& d) const
{
    write(OFstream(triSurfName(d))());
}


template<class Face>
void Foam::MeshedSurface<Face>::writeStats(Ostream& os) const
{
//    os  << "Faces        : " << size() << endl
//        << "Edges        : " << nEdges() << endl
//        << "Vertices     : " << nPoints() << endl
//        << "Bounding Box : " << boundBox(localPoints(), false) << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::operator=(const MeshedSurface& surf)
{
    clear();

    faces()  = surf.faces();
    points() = surf.points();
    patches_ = surf.patches_;
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MeshedSurfaceCleanup.C"
#include "MeshedSurfaceIO.C"
#include "MeshedSurfaceNew.C"

#undef  nativeSurfaceExt

// ************************************************************************* //

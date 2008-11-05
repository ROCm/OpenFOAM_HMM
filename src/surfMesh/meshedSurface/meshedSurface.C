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

#include "meshedSurface.H"
#include "keyedSurface.H"
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

namespace Foam
{

defineTypeNameAndDebug(meshedSurface, 0);

defineMemberFunctionSelectionTable
(
    meshedSurface,
    write,
    fileExtension
);

}

//  File extension for 'native' raw format
//! @cond localscope
const char * const nativeExt = "ofs";
//! @endcond localscope

// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

Foam::fileName Foam::meshedSurface::triSurfInstance(const Time& d)
{
    fileName foamName(d.caseName() + "." + nativeExt);

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
                    Pout<< " meshedSurface::triSurfInstance(const Time& d)"
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
        Pout<< " meshedSurface::triSurfInstance(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }

    return "constant";
}


Foam::fileName Foam::meshedSurface::triSurfName(const Time& d)
{
    fileName foamName(d.caseName() + "." + nativeExt);

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
                    Pout<< " meshedSurface::triSurfName(const Time& d)"
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
        Pout<< " meshedSurface::triSurfName(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }

    return d.path()/"constant"/typeName/foamName;
}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshedSurface::onePatch()
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


void Foam::meshedSurface::checkPatches()
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
                "meshedSurface::checkPatches()\n"
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
                "meshedSurface::checkPatches()\n"
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


void Foam::meshedSurface::sortFacesByRegion
(
    const List<label>& regionIds,
    const Map<word>& regionNames
)
{
    const List<FaceType>& unsortedFaces = faces();

    if (!&regionNames || !&regionIds || regionIds.size() == 0)
    {
        onePatch();
    }
    else if (regionIds.size() == unsortedFaces.size())
    {
        labelList faceMap;
        surfGroupList newPatches = keyedSurface::sortedRegions
        (
            regionIds,
            regionNames,
            faceMap
        );
        patches_.transfer(newPatches);

        // this is somewhat like ListOps reorder and/or IndirectList
        List<FaceType> newFaces(unsortedFaces.size());
        forAll(newFaces, faceI)
        {
            newFaces[faceI] = unsortedFaces[faceMap[faceI]];
        }
        faceMap.clear();

        faces().transfer(newFaces);
    }
}


// Read surf grouping, points, faces directly from Istream
bool Foam::meshedSurface::read(Istream& is, const bool doTriangulate)
{
    List<surfGroup> patchLst(is);
    is >> points() >> faces();

    // copy patch info:
    patches_.setSize(patchLst.size());
    forAll(patchLst, patchI)
    {
        patches_[patchI] = surfGroup
        (
            patchLst[patchI],
            patchI
        );
    }

    if (doTriangulate)
    {
        triangulate();
    }

    return is.good();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshedSurface::meshedSurface()
:
    MeshStorage(List<FaceType>(), pointField())
{}


Foam::meshedSurface::meshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<FaceType> >& faceLst,
    const xfer<surfGroupList>& patchLst
)
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_(patchLst)
{
    points().transfer(pointLst());
    faces().transfer(faceLst());
}


Foam::meshedSurface::meshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<FaceType> >& faceLst,
    const List<label>& patchSizes,
    const List<word>& patchNames,
    const List<word>& patchTypes
)
:
    MeshStorage(List<FaceType>(), pointField())
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


Foam::meshedSurface::meshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<FaceType> >& faceLst,
    const List<label>& regionIds,
    const Map<word>& regionNames
)
:
    MeshStorage(List<FaceType>(), pointField())
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
            "meshedSurface::meshedSurface(\n"
            "(\n"
            "    const xfer<pointField>&,\n"
            "    const xfer<List<FaceType> >&,\n"
            "    const List<label>& regionIds,\n"
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


Foam::meshedSurface::meshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<FaceType> >& faceLst,
    const List<label>& regionIds,
    const HashTable<label>& nameToRegionMapping
)
:
    MeshStorage(List<FaceType>(), pointField())
{
    points().transfer(pointLst());
    faces().transfer(faceLst());

    if (regionIds.size() != nFaces())
    {
        FatalErrorIn
        (
            "meshedSurface::meshedSurface(\n"
            "(\n"
            "    const xfer<pointField>&,\n"
            "    const xfer<List<FaceType> >&,\n"
            "    const List<label>& regionIds,\n"
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


Foam::meshedSurface::meshedSurface
(
    const polyBoundaryMesh& bMesh,
    const bool useGlobalPoints
)
:
    MeshStorage(List<FaceType>(), pointField())
{
    const polyMesh& mesh = bMesh.mesh();
    const polyPatchList& bPatches = bMesh;
    const label nIntFaces = mesh.nInternalFaces();

    // Get patch for all of outside
    primitivePatch allBoundary
    (
        SubList<FaceType>
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
Foam::meshedSurface::meshedSurface
(
    const surfMesh& sMesh
)
:
    MeshStorage(List<FaceType>(sMesh.faces()), sMesh.points())
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


Foam::meshedSurface::meshedSurface
(
    const keyedSurface& surf
)
:
    MeshStorage(List<FaceType>(), surf.points())
{
    labelList faceMap;
    surfGroupList patchLst = surf.sortedRegions(faceMap);
    patches_.transfer(patchLst);

    const List<FaceType>& origFaces = surf.faces();
    List<FaceType> newFaces(origFaces.size());

    // this is somewhat like ListOps reorder and/or IndirectList
    forAll(newFaces, faceI)
    {
        newFaces[faceI] = origFaces[faceMap[faceI]];
    }

    faces().transfer(newFaces);
}


Foam::meshedSurface::meshedSurface
(
    const fileName& fName,
    const word& ext,
    const bool triangulate
)
:
    MeshStorage(List<FaceType>(), pointField())
{
    read(fName, ext, triangulate);
}


Foam::meshedSurface::meshedSurface
(
    const fileName& fName,
    const bool triangulate
)
:
    MeshStorage(List<FaceType>(), pointField())
{
    read(fName, fName.ext(), triangulate);
}


Foam::meshedSurface::meshedSurface(Istream& is, const bool triangulate)
:
    MeshStorage(List<FaceType>(), pointField())
{
    read(is, triangulate);
}


Foam::meshedSurface::meshedSurface(const Time& d)
:
    MeshStorage(List<FaceType>(), pointField())
{
    read(IFstream(triSurfName(d))());
    // setDefaultPatches();
}


Foam::meshedSurface::meshedSurface(const meshedSurface& surf)
:
    MeshStorage(surf.faces(), surf.points()),
    patches_(surf.patches_)
{}


Foam::meshedSurface::meshedSurface(const xfer<keyedSurface>& surf)
:
    MeshStorage(List<FaceType>(), pointField())
{
    transfer(surf());
}


Foam::meshedSurface::meshedSurface(const xfer<meshedSurface>& surf)
:
    MeshStorage(List<FaceType>(), pointField())
{
    transfer(surf());
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshedSurface::~meshedSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//- Move points
void Foam::meshedSurface::movePoints(const pointField& newPoints)
{
    // Remove all geometry dependent data
    clearTopology();

    // Adapt for new point position
    MeshStorage::movePoints(newPoints);

    // Copy new points
    points() = newPoints;
}


//- scale points
void Foam::meshedSurface::scalePoints(const scalar& scaleFactor)
{
    // avoid bad scaling
    if (scaleFactor > 0 && scaleFactor != 1.0)
    {
        // Remove all geometry dependent data
        clearTopology();

        // Adapt for new point position
        MeshStorage::movePoints(pointField());

        points() *= scaleFactor;
    }
}


Foam::meshedSurface Foam::meshedSurface::subsetMesh
(
    const boolList& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const pointField& locPoints = localPoints();
    const List<FaceType>& locFaces = localFaces();

    // Fill pointMap, faceMap
    subsetMap(include, pointMap, faceMap);

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
    List<FaceType> newFaces(faceMap.size());

    forAll(faceMap, faceI)
    {
        // Get old vertex labels
        label origFaceI = faceMap[faceI];
        const FaceType& oldFace = locFaces[origFaceI];

        newFaces[faceI] = FaceType(oldFace);

        // Renumber labels for face
        FaceType& f = newFaces[faceI];
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

    // adjust patch start
    label startFaceI = 0;
    forAll(newPatches, patchI)
    {
        newPatches[patchI].start() = startFaceI;
        startFaceI += newPatches[patchI].size();
    }

    // Construct an empty subsurface and fill
    meshedSurface subSurf;

    subSurf.patches_.transfer(newPatches);
    subSurf.points().transfer(newPoints);
    subSurf.faces().transfer(newFaces);

    return subSurf;
}


void Foam::meshedSurface::transfer(meshedSurface& surf)
{
    clearOut();
    points().transfer(surf.points());
    faces().transfer(surf.faces());
    patches_.transfer(surf.patches_);
    surf.clearOut();
}


void Foam::meshedSurface::transfer(keyedSurface& surf)
{
    clearOut();
    points().transfer(surf.points());
    faces().clear();

    labelList faceMap;
    surfGroupList patchLst = surf.sortedRegions(faceMap);
    patches_.transfer(patchLst);
    surf.regions().clear();
    surf.patches_.clear();

    List<FaceType>& oldFaces = surf.faces();
    List<FaceType> newFaces(oldFaces.size());

    // this is somewhat like ListOps reorder and/or IndirectList
    forAll(newFaces, faceI)
    {
        newFaces[faceI] = oldFaces[faceMap[faceI]];
    }

    faces().transfer(newFaces);
    surf.faces().clear();
    surf.clearOut();
}


bool Foam::meshedSurface::canRead(const word& ext, const bool verbose)
{
    return keyedSurface::canRead(ext, verbose);
}


bool Foam::meshedSurface::canWrite(const word& ext, const bool verbose)
{
    // perhaps we got sent an entire name
    word fExt(ext);

    // FIXME: this looks horrible, but I don't have STL docs here
    string::size_type dot = ext.find_last_of(".");
    if (dot != string::npos)
    {
        fExt = ext.substr(dot+1);
    }

    // handle 'native' format directly
    if (fExt == nativeExt)
    {
        return true;
    }

    writefileExtensionMemberFunctionTable::iterator mfIter =
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
            Info<<"Valid types: ( " << nativeExt;
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


// Read from file in given format
bool Foam::meshedSurface::read
(
    const fileName& fName,
    const word& ext,
    const bool triangulate
)
{
    // handle 'native' format directly
    if (ext == nativeExt)
    {
        return read(IFstream(fName)(), triangulate);
    }
    else
    {
        // use selector mechanism
        transfer(New(fName, ext, triangulate)());
        return true;
    }
}


void Foam::meshedSurface::write
(
    const fileName& fName,
    const meshedSurface& surf
)
{
    if (debug)
    {
        Info<< "meshedSurface::write(const fileName&, const meshedSurface&) : "
               "writing meshedSurface to " << fName
            << endl;
    }

    const word ext = fName.ext();

    // handle 'native' format directly
    if (ext == nativeExt)
    {
        surf.write(OFstream(fName)());
    }
    else
    {
        writefileExtensionMemberFunctionTable::iterator mfIter =
            writefileExtensionMemberFunctionTablePtr_->find(ext);

        if (mfIter == writefileExtensionMemberFunctionTablePtr_->end())
        {
            FatalErrorIn
            (
                "meshedSurface::write(const fileName&)"
            )   << "Unknown file extension " << ext << nl << nl
                << "Valid types are :" << endl
                << writefileExtensionMemberFunctionTablePtr_->toc()
                << exit(FatalError);
        }

        mfIter()(fName, surf);
    }
}


void Foam::meshedSurface::write(Ostream& os) const
{
    // just emit some information until we get a nice IOobject
    IOobject::writeBanner(os);
    os  << "// OpenFOAM Surface format" << nl
        << "// ~~~~~~~~~~~~~~~~~~~~~~~" << nl
        << "// regions:" << nl
        << patches_.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    forAll(patches_, patchI)
    {
        patches_[patchI].writeDict(os);
    }
    os  << decrIndent << token::END_LIST << nl;

    IOobject::writeDivider(os);

    // Note: Write with global point numbering
    os  << "\n// points:" << nl << points() << nl;

    IOobject::writeDivider(os);
    os  << "\n// faces:"  << nl << faces() << nl;

    IOobject::writeDivider(os);

    // Check state of Ostream
    os.check("meshedSurface::write(Ostream&)");
}


void Foam::meshedSurface::write(const Time& d) const
{
    write(OFstream(triSurfName(d))());
}


void Foam::meshedSurface::writeStats(Ostream& os) const
{
    os  << "Faces        : " << size() << endl
        << "Edges        : " << nEdges() << endl
        << "Vertices     : " << nPoints() << endl
        << "Bounding Box : " << boundBox(localPoints(), false) << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::meshedSurface::operator=(const meshedSurface& surf)
{
    clearOut();
    faces()  = surf.faces();
    points() = surf.points();
    patches_ = surf.patches_;
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const meshedSurface& surf)
{
    surf.write(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

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

Foam::word Foam::meshedSurface::defaultGeometricType("empty");
//  File extension for 'native' raw format
//! @cond localscope
const char * const nativeExt = "ftr";
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
    patches_[0] = surfacePatch
    (
        defaultGeometricType,
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


// Read triangles, points from Istream
bool Foam::meshedSurface::read(Istream& is)
{
    is  >> patches_ >> points() >> faces();

    return true;
}

#if 0
// Read from file in given format
bool Foam::meshedSurface::read(const fileName& name, const word& ext)
{
    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();

        return read(unzipName, unzipName.ext());
    }
}
#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshedSurface::meshedSurface()
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
{}


Foam::meshedSurface::meshedSurface
(
    const pointField& pointLst,
    const List<FaceType>& faceLst,
    const surfacePatchList& patchLst
)
:
    MeshStorage(faceLst, pointLst),
    patches_(patchLst)
{}


Foam::meshedSurface::meshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<FaceType> >& faceLst,
    const xfer<surfacePatchList>& patchLst
)
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_(patchLst)
{
    faces().transfer(faceLst());
    points().transfer(pointLst());
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
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
{
    points().transfer(pointLst());
    faces().transfer(faceLst());

    surfacePatchList newPatches(patchSizes.size());

    label start = 0;
    forAll(newPatches, patchI)
    {
        newPatches[patchI] = surfacePatch
        (
            defaultGeometricType,
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
    const List<label>& regionLst,
    const Map<word>& regionNames
)
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
{
    points().transfer(pointLst());
    faces().transfer(faceLst());

    const List<FaceType>& unsortedFaces = faces();

    if (regionLst.size() == 0)
    {
        onePatch();
    }
    else if (regionLst.size() != unsortedFaces.size())
    {
        FatalErrorIn
        (
            "meshedSurface::meshedSurface(\n"
            "(\n"
            "    const pointField&,\n"
            "    const List<FaceType>&,\n"
            "    const List<label>& regionLst,\n"
            "    const Map<word>& regionNames\n"
            " )\n"
        )
            << "size mismatch : regionLst.size() != faceLst.size()"
            << exit(FatalError);
    }
    else
    {
        labelList faceMap;
        surfacePatchList newPatches = keyedSurface::sortedRegions
        (
            regionLst,
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


Foam::meshedSurface::meshedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<FaceType> >& faceLst
)
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
{
    faces().transfer(faceLst());
    points().transfer(pointLst());

    onePatch();
}


Foam::meshedSurface::meshedSurface
(
    const xfer<pointField>& pointLst,
    const List<keyedFace>& faceLst,
    const Map<word>& regionNames
)
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
{
    points().transfer(pointLst());

    labelList faceMap;
    surfacePatchList newPatches = keyedSurface::sortedRegions
    (
        faceLst,
        regionNames,
        faceMap
    );
    patches_.transfer(newPatches);

    // this is somewhat like ListOps reorder and/or IndirectList
    List<FaceType> newFaces(faceLst.size());
    forAll(newFaces, faceI)
    {
        newFaces[faceI] = faceLst[faceMap[faceI]];
    }

    faces().transfer(newFaces);
}


Foam::meshedSurface::meshedSurface
(
    const xfer<pointField>& pointLst,
    const List<keyedFace>& faceLst,
    const HashTable<label>& nameToRegionMapping
)
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
{
    points().transfer(pointLst());

    Map<word> regionNames;
    forAllConstIter(HashTable<label>, nameToRegionMapping, iter)
    {
        regionNames.insert(iter(), iter.key());
    }


    labelList faceMap;
    surfacePatchList newPatches = keyedSurface::sortedRegions
    (
        faceLst,
        regionNames,
        faceMap
    );
    patches_.transfer(newPatches);

    List<FaceType> newFaces(faceLst.size());

    // this is somewhat like ListOps reorder and/or IndirectList
    forAll(newFaces, faceI)
    {
        newFaces[faceI] = faceLst[faceMap[faceI]];
    }

    faces().transfer(newFaces);
}


Foam::meshedSurface::meshedSurface
(
    const polyBoundaryMesh& bMesh,
    const bool useGlobalPoints
)
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
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
    surfacePatchList newPatches(bPatches.size());

    label startFaceI = 0;
    forAll(bPatches, patchI)
    {
        const polyPatch& p = bPatches[patchI];

        newPatches[patchI] = surfacePatch
        (
            defaultGeometricType,
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
    MeshStorage(List<FaceType>(sMesh.faces()), sMesh.points()),
    patches_()
{
    const surfPatchList& sPatches = sMesh.boundaryMesh();

    // create patch list
    surfacePatchList newPatches(sPatches.size());

    label startFaceI = 0;
    forAll(sPatches, patchI)
    {
        const surfPatch& p = sPatches[patchI];

        newPatches[patchI] = surfacePatch
        (
            defaultGeometricType,
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
    MeshStorage(List<FaceType>(), surf.points()),
    patches_()
{
    labelList faceMap;
    surfacePatchList newPatches = surf.sortedRegions(faceMap);
    patches_.transfer(newPatches);

    const List<keyedFace>& origFaces = surf.faces();
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
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
{
    // use selector mechanism
    autoPtr<meshedSurface> surfPtr = New(fName, ext, triangulate);
    transfer(surfPtr());
}

Foam::meshedSurface::meshedSurface
(
    const fileName& fName,
    const bool triangulate
)
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
{
    // use selector mechanism
    autoPtr<meshedSurface> surfPtr = New(fName, triangulate);
    transfer(surfPtr());
}


Foam::meshedSurface::meshedSurface(Istream& is)
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
{
    read(is);
    // setDefaultPatches();
}


Foam::meshedSurface::meshedSurface(const Time& d)
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
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
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
{
    transfer(surf());
}


Foam::meshedSurface::meshedSurface(const xfer<meshedSurface>& surf)
:
    MeshStorage(List<FaceType>(), pointField()),
    patches_()
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
    surfacePatchList newPatches(patches_);
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
    faces().transfer(surf.faces());
    points().transfer(surf.points());
    patches_.transfer(surf.patches_);
    surf.clearOut();
}


void Foam::meshedSurface::transfer(keyedSurface& surf)
{
    clearOut();
    points().transfer(surf.points());

    labelList faceMap;
    surfacePatchList newPatches = surf.sortedRegions(faceMap);
    patches_.transfer(newPatches);

    const List<keyedFace>& origFaces = surf.faces();
    List<FaceType> newFaces(origFaces.size());

    // this is somewhat like ListOps reorder and/or IndirectList
    forAll(newFaces, faceI)
    {
        newFaces[faceI] = origFaces[faceMap[faceI]];
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
            const wordList& known =
                writefileExtensionMemberFunctionTablePtr_->toc();

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
    // quick-hack
    os << patches_.size() << nl << token::BEGIN_LIST;
    forAll(patches_, patchI)
    {
        patches_[patchI].writeDict(os);
    }
    os << token::END_LIST;

    // Note: Write with global point numbering
    os << points() << nl << faces() << endl;

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

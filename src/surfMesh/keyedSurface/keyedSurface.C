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

#include "keyedSurface.H"
#include "meshedSurface.H"
#include "demandDrivenData.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
#include "boundBox.H"
#include "SortableList.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(keyedSurface, 0);
defineRunTimeSelectionTable(keyedSurface, fileExtension);
defineMemberFunctionSelectionTable
(
    keyedSurface,
    write,
    fileExtension
);
}

Foam::word Foam::keyedSurface::defaultGeometricType("empty");

//  File extension for 'native' raw format
//! @cond localscope
const char * const nativeExt = "ftr";
//! @endcond localscope

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::keyedSurface::triSurfInstance(const Time& d)
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
                    Pout<< " keyedSurface::triSurfInstance(const Time& d)"
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
        Pout<< " keyedSurface::triSurfInstance(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }

    return "constant";
}


Foam::fileName Foam::keyedSurface::triSurfName(const Time& d)
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
                    Pout<< " keyedSurface::triSurfName(const Time& d)"
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
        Pout<< " keyedSurface::triSurfName(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }

    return d.path()/"constant"/typeName/foamName;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::string Foam::keyedSurface::getLineNoComment(IFstream& is)
{
    string line;
    do
    {
        is.getLine(line);
    }
    while ((line.size() == 0 || line[0] == '#') && is.good());

    return line;
}


void Foam::keyedSurface::setPatches(const label maxPatch)
{
    geoPatches_.setSize(maxPatch+1);

    forAll(geoPatches_, patchI)
    {
        geoPatches_[patchI] = geometricSurfacePatch
        (
            defaultGeometricType,
            "patch" + ::Foam::name(patchI),
            patchI
        );
    }
}


void Foam::keyedSurface::setPatches()
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



void Foam::keyedSurface::setPatches
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

    geoPatches_.setSize(maxPatch+1);

    forAll(geoPatches_, patchI)
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

        geoPatches_[patchI] = geometricSurfacePatch
        (
            defaultGeometricType,
            patchName,
            patchI
        );
    }
}


void Foam::keyedSurface::setPatches
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

// Read triangles, points from Istream
//
//  FIXME: needs to handle labelledTri etc.
bool Foam::keyedSurface::read(Istream& is)
{
    notImplemented("Foam::keyedSurface::read(Istream&)");
    return false;

    //// is  >> geoPatches_ >> points() >> faces();
    //// return true;
}

#if 0
// Read from file in given format
bool Foam::keyedSurface::read(const fileName& name, const word& ext)
{
    if (ext == "gz")
    {
        fileName unzipName = name.lessExt();

        return read(unzipName, unzipName.ext());
    }
}
#endif


// Returns patch info.
// Sets faceMap to the indexing according to patch numbers.
// Patch numbers start at 0.
Foam::surfacePatchList Foam::keyedSurface::sortedRegions
(
    const List<label>& regionLst,
    const Map<word>& patchNames,
    labelList& faceMap
)
{
    // determine sort order according to region numbers

    // std::sort() really seems to mix up the order.
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
    surfacePatchList surfPatches(regionLookup.size());
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

        surfPatches[patchI] = surfacePatch
        (
            defaultGeometricType,
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

        faceMap[faceI] =
            surfPatches[patchI].start() + surfPatches[patchI].size()++;
    }

    // with reordered faces
    return surfPatches;
}


Foam::surfacePatchList Foam::keyedSurface::sortedRegions
(
    labelList& faceMap
) const
{
    // supply some patch names
    Map<word> patchNames;
    forAll(geoPatches_, patchI)
    {
        patchNames.insert(patchI, geoPatches_[patchI].name());
    }

    return sortedRegions(regions_, patchNames, faceMap);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::keyedSurface::keyedSurface()
:
    MeshStorage(List<FaceType>(), pointField())
{}


Foam::keyedSurface::keyedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<face> >& faceLst,
    const xfer<List<label> >& regionIds,
    const xfer<geometricSurfacePatchList>& patchLst
)
:
    MeshStorage(List<FaceType>(), pointField()),
    regions_(regionIds),
    geoPatches_(patchLst)
{
    points().transfer(pointLst());
    faces().transfer(faceLst());
}


Foam::keyedSurface::keyedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<face> >& faceLst,
    const xfer<List<label> >& regionIds,
    const Map<word>& regionNames
)
:
    MeshStorage(List<FaceType>(), pointField()),
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


Foam::keyedSurface::keyedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<face> >& faceLst,
    const xfer<List<label> >& regionIds,
    const HashTable<label>& labelToRegion
)
:
    MeshStorage(List<FaceType>(), pointField()),
    regions_(regionIds)
{
    points().transfer(pointLst());
    faces().transfer(faceLst());

    // set patch names from (name => id) mapping
    setPatches(labelToRegion);
}


Foam::keyedSurface::keyedSurface
(
    const xfer<pointField>& pointLst,
    const xfer<List<face> >& faceLst
)
:
    MeshStorage(List<FaceType>(), pointField()),
    regions_(faceLst().size(), 0),     // single default patch
    geoPatches_()
{
    points().transfer(pointLst());
    faces().transfer(faceLst());

    setPatches(0);
}



Foam::keyedSurface::keyedSurface
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

    geometricSurfacePatchList newPatches(bPatches.size());

    // Get patch for all of outside
    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nFaces() - nIntFaces,
            nIntFaces
        ),
        mesh.points()
    );

    List<FaceType> newFaces(allBoundary.size());
    List<label>    newRegions(allBoundary.size());

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
    const List<face>& bfaces =
    (
        useGlobalPoints
      ? allBoundary
      : allBoundary.localFaces()
    );

    label faceIndex = 0;
    forAll(bPatches, patchI)
    {
        const polyPatch& p = bPatches[patchI];

        newPatches[patchI] = geometricSurfacePatch
        (
            defaultGeometricType,
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
    geoPatches_.transfer(newPatches);
}


Foam::keyedSurface::keyedSurface
(
    const meshedSurface& surf
)
:
    MeshStorage(List<FaceType>(), surf.points())
{
    const List<face>& origFaces = surf.faces();
    const surfacePatchList& patchLst = surf.patches();

    List<FaceType> newFaces(origFaces.size());
    List<label>    newRegions(origFaces.size());
    geometricSurfacePatchList newPatches(patchLst.size());

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
    geoPatches_.transfer(newPatches);
}


Foam::keyedSurface::keyedSurface
(
    const fileName& fName,
    const word& ext,
    const bool triangulate
)
:
    MeshStorage(List<FaceType>(), pointField())
{
    // use selector mechanism
    autoPtr<keyedSurface> surfPtr = New(fName, ext, triangulate);
    transfer(surfPtr());
}


Foam::keyedSurface::keyedSurface
(
    const fileName& fName,
    const bool triangulate
)
:
    MeshStorage(List<FaceType>(), pointField())
{
    // use selector mechanism
    autoPtr<keyedSurface> surfPtr = New(fName, triangulate);
    transfer(surfPtr());
}


Foam::keyedSurface::keyedSurface(Istream& is)
:
    MeshStorage(List<FaceType>(), pointField())
{
    read(is);
    setPatches();
}


Foam::keyedSurface::keyedSurface(const Time& d)
:
    MeshStorage(List<FaceType>(), pointField())
{
    read(IFstream(triSurfName(d))());
    setPatches();
}


Foam::keyedSurface::keyedSurface(const keyedSurface& surf)
:
    MeshStorage(surf.faces(), surf.points()),
    regions_(surf.regions_),
    geoPatches_(surf.geoPatches_)
{}


Foam::keyedSurface::keyedSurface(const xfer<keyedSurface>& surf)
:
    MeshStorage(List<FaceType>(), pointField()),
    regions_(),
    geoPatches_()
{
    transfer(surf());
}


Foam::keyedSurface::keyedSurface(const xfer<meshedSurface>& surf)
:
    MeshStorage(List<FaceType>(), pointField()),
    regions_(),
    geoPatches_()
{
    transfer(surf());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::keyedSurface::~keyedSurface()
{
    // Info<<"~keyedSurface()" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Move points
void Foam::keyedSurface::movePoints(const pointField& newPoints)
{
    // Remove all geometry dependent data
    clearTopology();

    // Adapt for new point position
    MeshStorage::movePoints(newPoints);

    // Copy new points
    points() = newPoints;
}


//- scale points
void Foam::keyedSurface::scalePoints(const scalar& scaleFactor)
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



Foam::keyedSurface Foam::keyedSurface::subsetMesh
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

    // Renumber face node labels and compact
    List<FaceType> newFaces(faceMap.size());

    forAll(faceMap, faceI)
    {
        // Get old vertex labels
        const FaceType& oldFace = locFaces[faceMap[faceI]];

        newFaces[faceI] = FaceType(oldFace);

        // Renumber labels for face
        FaceType& f = newFaces[faceI];
        forAll(f, fp)
        {
            f[fp] = oldToNew[oldFace[fp]];
        }
    }

    // Construct an empty subsurface and fill
    keyedSurface subSurf;

    // transfer
    subSurf.points().transfer(newPoints);
    subSurf.faces().transfer(newFaces);

    // copy
    subSurf.geoPatches_ = geoPatches_;

    return subSurf;
}


void Foam::keyedSurface::transfer(keyedSurface& surf)
{
    clearOut();
    faces().transfer(surf.faces());
    points().transfer(surf.points());
    regions_.transfer(surf.regions_);
    geoPatches_.transfer(surf.geoPatches_);
    surf.clearOut();
}


void Foam::keyedSurface::transfer(meshedSurface& surf)
{
    clearOut();
    points().transfer(surf.points());
    faces().transfer(surf.faces());
    regions_.setSize(nFaces());

    surfacePatchList& patchLst = surf.patches();

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        // copy info
        geoPatches_[patchI] = patchLst[patchI];

        forAll(patchLst[patchI], patchFaceI)
        {
            regions_[faceIndex++] = patchI;
        }
    }

    patchLst.clear();
    surf.clearOut();
}


bool Foam::keyedSurface::canRead(const word& ext, const bool verbose)
{
    // FIXME: this looks horrible, but I don't have my STL docs here

    // perhaps we got sent an entire name
    word fExt(ext);

    string::size_type dot = fExt.find_last_of(".");
    if (dot != string::npos)
    {
        fExt = fExt.substr(dot+1);
    }

    fileExtensionConstructorTable::iterator cstrIter =
        fileExtensionConstructorTablePtr_->find(fExt);

    // handle 'native' format directly
    if (fExt == nativeExt)
    {
        return true;
    }

    // would be nice to have information about which format this actually is
    if (cstrIter == fileExtensionConstructorTablePtr_->end())
    {
        if (verbose)
        {
            const wordList& known = fileExtensionConstructorTablePtr_->toc();

            Info<<"Unknown file extension for reading: " << fExt << nl;
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


bool Foam::keyedSurface::canWrite(const word& ext, const bool verbose)
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



#if 0
void Foam::keyedSurface::write
(
    const fileName& fName,
    const bool sorted
) const
{
    write(fName, fName.ext(), sorted);
}
#endif


void Foam::keyedSurface::write
(
    const fileName& fName,
    const keyedSurface& surf
)
{
    if (debug)
    {
        Info<< "keyedSurface::write(const fileName&, const keyedSurface&) : "
               "writing keyedSurface to " << fName
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
                "keyedSurface::write(const fileName&)"
            )   << "Unknown file extension " << ext << nl << nl
                << "Valid types are :" << endl
                << writefileExtensionMemberFunctionTablePtr_->toc()
                << exit(FatalError);
        }

        mfIter()(fName, surf);
    }
}



void Foam::keyedSurface::write(Ostream& os) const
{
    os  << "\n// geometric regions:\n"
        << geoPatches_ << endl;

    // Note: Write with global point numbering
    os  << "\n// points:"
        << points() << nl
        << "\n// faces:"
        << faces() << nl
        << "\n// regions:"
        << regions() << endl;

    // Check state of Ostream
    os.check("keyedSurface::write(Ostream&)");
}


void Foam::keyedSurface::write(const Time& d) const
{
    write(OFstream(triSurfName(d))());
}


void Foam::keyedSurface::writeStats(Ostream& os) const
{
    os  << "Faces        : " << size() << endl
        << "Edges        : " << nEdges() << endl
        << "Vertices     : " << nPoints() << endl
        << "Bounding Box : " << boundBox(localPoints(), false) << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::keyedSurface::operator=(const keyedSurface& surf)
{
    clearOut();
    faces()  = surf.faces();
    points() = surf.points();
    geoPatches_ = surf.geoPatches_;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const keyedSurface& surf)
{
    surf.write(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

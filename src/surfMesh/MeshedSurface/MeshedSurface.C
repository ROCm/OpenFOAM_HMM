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
#include "surfMesh.H"
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
    const Xfer<surfZoneList>& zoneLst
)
:
    ParentType(pointLst, faceLst),
    zones_(zoneLst)
{}


template<class Face>
Foam::MeshedSurface<Face>::MeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst,
    const UList<label>& zoneSizes,
    const UList<word>& zoneNames
)
:
    ParentType(pointLst, faceLst)
{
    if (&zoneSizes)
    {
        if (&zoneNames)
        {
            addZones(zoneSizes, zoneNames);
        }
        else
        {
            addZones(zoneSizes);
        }
    }
    else
    {
        oneZone();
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

    label startFaceI = 0;
    label nZone = 0;
    forAll(bPatches, patchI)
    {
        const polyPatch& p = bPatches[patchI];

        if (p.size())
        {
            newZones[nZone] = surfZone
            (
                p.name(),
                p.size(),
                startFaceI,
                nZone
            );

            nZone++;
            startFaceI += p.size();
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
Foam::MeshedSurface<Face>::MeshedSurface(const surfMesh& mesh)
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
    const UnsortedMeshedSurface<Face>& surf
)
{
    labelList faceMap;
    surfZoneList zoneLst = surf.sortedZones(faceMap);
    zones_.transfer(zoneLst);

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
    zones_(surf.zones_)
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::oneZone(const word& name)
{
    word zoneName(name);
    if (zoneName.empty())
    {
        if (zones_.size())
        {
            zoneName = zones_[0].name();
        }
        if (zoneName.empty())
        {
            zoneName = "zone0";
        }
    }

    // set single default zone
    zones_.setSize(1);
    zones_[0] = surfZone
    (
        zoneName,
        this->size(),   // zone size
        0,              // zone start
        0               // zone index
    );
}


template<class Face>
void Foam::MeshedSurface<Face>::checkZones()
{
    // extra safety, ensure we have at some zones
    // and they cover all the faces - fix start silently
    if (zones_.size() <= 1)
    {
        oneZone();
    }
    else
    {
        label count = 0;
        forAll(zones_, zoneI)
        {
            zones_[zoneI].start() = count;
            count += zones_[zoneI].size();
        }

        if (count < size())
        {
            WarningIn
            (
                "MeshedSurface::checkZones()\n"
            )
                << "more faces " << size() << " than zones " << count
                << " ... extending final zone"
                << endl;

            zones_[zones_.size()-1].size() += count - size();
        }
        else if (count > size())
        {
            FatalErrorIn
            (
                "MeshedSurface::checkZones()\n"
            )
                << "more zones " << count << " than faces " << size()
                << exit(FatalError);
        }
    }
}


template<class Face>
void Foam::MeshedSurface<Face>::sortFacesAndStore
(
    const Xfer<List<Face> >& unsortedFaces,
    const Xfer<List<label> >& zoneIds,
    const bool sorted
)
{
    List<Face>  oldFaces(unsortedFaces);
    List<label> zones(zoneIds);

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
        sortedOrder(zones, faceMap);
        zones.clear();

        // sorted faces
        List<Face> newFaces(faceMap.size());
        forAll(faceMap, faceI)
        {
            // use transfer to recover memory if possible
            newFaces[faceI].transfer(oldFaces[faceMap[faceI]]);
        }
        this->storedFaces().transfer(newFaces);

    }
    zones.clear();
}


template<class Face>
void Foam::MeshedSurface<Face>::remapFaces
(
    const UList<label>& faceMap
)
{
    // recalculate the zone start/size
    if (&faceMap && faceMap.size())
    {
        if (zones_.empty())
        {
            oneZone();
        }
        else if (zones_.size() == 1)
        {
            // optimized for single zone case
            zones_[0].size() = faceMap.size();
        }
        else
        {
            label newFaceI = 0;
            label origEndI = 0;
            forAll(zones_, zoneI)
            {
                surfZone& zone = zones_[zoneI];

                // adjust zone start
                zone.start() = newFaceI;
                origEndI += zone.size();

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

                // adjust zone size
                zone.size() = newFaceI - zone.start();
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::clear()
{
    ParentType::clear();
    zones_.clear();
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

    // create/copy a new zones list, each zone with zero size
    surfZoneList newZones(zones_);
    forAll(newZones, zoneI)
    {
        newZones[zoneI].size() = 0;
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

    // recalculate the zones start/size
    label newFaceI = 0;
    label origEndI = 0;

    // adjust zone sizes
    forAll(newZones, zoneI)
    {
        surfZone& zone = newZones[zoneI];

        // adjust zone start
        zone.start() = newFaceI;
        origEndI += zone.size();

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

        // adjust zone size
        zone.size() = newFaceI - zone.start();
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
void Foam::MeshedSurface<Face>::addZones
(
    const UList<surfZone>& zones,
    const bool cullEmpty
)
{
    label nZone = 0;

    zones_.setSize(zones.size());
    forAll(zones_, zoneI)
    {
        if (zones[zoneI].size() || !cullEmpty)
        {
            zones_[nZone] = surfZone(zones[zoneI], nZone);
            nZone++;
        }
    }
    zones_.setSize(nZone);
}


template<class Face>
void Foam::MeshedSurface<Face>::addZones
(
    const UList<label>& sizes,
    const UList<word>& names,
    const bool cullEmpty
)
{
    label start   = 0;
    label nZone = 0;

    zones_.setSize(sizes.size());
    forAll(zones_, zoneI)
    {
        if (sizes[zoneI] || !cullEmpty)
        {
            zones_[nZone] = surfZone
            (
                names[zoneI],
                sizes[zoneI],
                start,
                nZone
            );
            start += sizes[zoneI];
            nZone++;
        }
    }
    zones_.setSize(nZone);
}


template<class Face>
void Foam::MeshedSurface<Face>::addZones
(
    const UList<label>& sizes,
    const bool cullEmpty
)
{
    label start   = 0;
    label nZone = 0;

    zones_.setSize(sizes.size());
    forAll(zones_, zoneI)
    {
        if (sizes[zoneI] || !cullEmpty)
        {
            zones_[nZone] = surfZone
            (
                word("zone") + ::Foam::name(nZone),
                sizes[zoneI],
                start,
                nZone
            );
            start += sizes[zoneI];
            nZone++;
        }
    }
    zones_.setSize(nZone);
}


template<class Face>
void Foam::MeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    reset(xferMove(surf.storedPoints()), xferMove(surf.storedFaces()));
    zones_.transfer(surf.zones_);

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
    List<Face>& oldFaces = surf.storedFaces();

    List<Face> newFaces(faceMap.size());
    forAll(faceMap, faceI)
    {
        newFaces[faceI].transfer(oldFaces[faceMap[faceI]]);
    }
    faceMap.clear();

    reset(xferMove(surf.storedPoints()), xferMove(newFaces));
    zones_.transfer(zoneLst);

    surf.clear();
}


template<class Face>
Foam::Xfer< Foam::MeshedSurface<Face> >
Foam::MeshedSurface<Face>::xfer()
{
    return xferMove(*this);
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
    zones_ = surf.zones_;
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MeshedSurfaceIO.C"
#include "MeshedSurfaceNew.C"

// ************************************************************************* //

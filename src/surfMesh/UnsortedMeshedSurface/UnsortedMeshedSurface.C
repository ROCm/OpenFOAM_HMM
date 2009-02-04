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
#include "polyBoundaryMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
Foam::wordHashSet Foam::UnsortedMeshedSurface<Face>::readTypes()
{
    return wordHashSet(*fileExtensionConstructorTablePtr_);
}


template<class Face>
Foam::wordHashSet Foam::UnsortedMeshedSurface<Face>::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
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
void Foam::UnsortedMeshedSurface<Face>::write
(
    const fileName& name,
    const UnsortedMeshedSurface<Face>& surf
)
{
    if (debug)
    {
        Info<< "UnsortedMeshedSurface::write"
            "(const fileName&, const UnsortedMeshedSurface&) : "
            "writing to " << name
            << endl;
    }

    const word ext = name.ext();

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
            "UnsortedMeshedSurface::write"
            "(const fileName&, const UnsortedMeshedSurface&)"
        )   << "Unknown file extension " << ext << nl << nl
            << "Valid types are :" << endl
            << writeTypes()
            << exit(FatalError);
    }

    mfIter()(name, surf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface()
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst,
    const Xfer<List<label> >& zoneIds,
    const Xfer<surfZoneIdentifierList>& zoneTofc
)
:
    ParentType(pointLst, faceLst),
    zoneIds_(zoneIds),
    zoneToc_(zoneTofc)
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
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
            setZones(zoneSizes, zoneNames);
        }
        else
        {
            setZones(zoneSizes);
        }
    }
    else
    {
        oneZone();
    }
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
    setZones(surf.zones());
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const fileName& name,
    const word& ext
)
{
    read(name, ext);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface(const fileName& name)
{
    read(name);
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
    zoneIds_(surf.zoneIds_),
    zoneToc_(surf.zoneToc_)
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const Xfer<UnsortedMeshedSurface<Face> >& surf
)
{
    transfer(surf());
}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const Xfer<MeshedSurface<Face> >& surf
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
void Foam::UnsortedMeshedSurface<Face>::oneZone(const word& name)
{
    zoneIds_.setSize(size());
    zoneIds_ = 0;

    word zoneName(name);
    if (zoneName.empty())
    {
        if (zoneToc_.size())
        {
            zoneName = zoneToc_[0].name();
        }
        if (zoneName.empty())
        {
            zoneName = "zone0";
        }
    }

    // set single default zone
    zoneToc_.setSize(1);
    zoneToc_[0] = surfZoneIdentifier(zoneName, 0);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setZones
(
    const surfZoneList& zoneLst
)
{
    zoneIds_.setSize(size());
    zoneToc_.setSize(zoneLst.size());

    forAll(zoneToc_, zoneI)
    {
        const surfZone& zone = zoneLst[zoneI];

        zoneToc_[zoneI] = zone;

        // assign sub-zone Ids
        SubList<label> subZone(zoneIds_, zone.size(), zone.start());
        subZone = zoneI;
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setZones
(
    const UList<label>& sizes,
    const UList<word>& names
)
{
    zoneIds_.setSize(size());
    zoneToc_.setSize(sizes.size());

    label start = 0;
    forAll(zoneToc_, zoneI)
    {
        zoneToc_[zoneI] = surfZoneIdentifier(names[zoneI], zoneI);

        // assign sub-zone Ids
        SubList<label> subZone(zoneIds_, sizes[zoneI], start);
        subZone = zoneI;

        start += sizes[zoneI];
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setZones
(
    const UList<label>& sizes
)
{
    zoneIds_.setSize(size());
    zoneToc_.setSize(sizes.size());

    label start = 0;
    forAll(zoneToc_, zoneI)
    {
        zoneToc_[zoneI] = surfZoneIdentifier
        (
            word("zone") + ::Foam::name(zoneI),
            zoneI
        );

        // assign sub-zone Ids
        SubList<label> subZone(zoneIds_, sizes[zoneI], start);
        subZone = zoneI;

        start += sizes[zoneI];
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::remapFaces
(
    const UList<label>& faceMap
)
{
    // re-assign the zone Ids
    if (&faceMap && faceMap.size())
    {
        if (zoneToc_.empty())
        {
            oneZone();
        }
        else if (zoneToc_.size() == 1)
        {
            // optimized for single-zone case
            zoneIds_ = 0;
        }
        else
        {
            List<label> newZones(faceMap.size());

            forAll(faceMap, faceI)
            {
                newZones[faceI] = zoneIds_[faceMap[faceI]];
            }
            zoneIds_.transfer(newZones);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setSize(const label s)
{
    ParentType::setSize(s);
    // if zones extend: set with last zoneId
    zoneIds_.setSize(s, zoneToc_.size() - 1);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::clear()
{
    ParentType::clear();
    zoneIds_.clear();
    zoneToc_.clear();
}


template<class Face>
Foam::surfZoneList Foam::UnsortedMeshedSurface<Face>::sortedZones
(
    labelList& faceMap
) const
{
    // supply some zone names
    Map<word> zoneNames;
    forAll(zoneToc_, zoneI)
    {
        zoneNames.insert(zoneI, zoneToc_[zoneI].name());
    }

    return sortedZonesById(zoneIds_, zoneNames, faceMap);
}


template<class Face>
Foam::UnsortedMeshedSurface<Face> Foam::UnsortedMeshedSurface<Face>::subsetMesh
(
    const labelHashSet& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const pointField&  locPoints = this->localPoints();
    const List<Face>&  locFaces  = this->localFaces();

    // Fill pointMap, faceMap
    PatchTools::subsetMap(*this, include, pointMap, faceMap);

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
    List<label> newZones(faceMap.size());

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

        newZones[faceI] = zoneIds_[origFaceI];
    }
    oldToNew.clear();

    // construct a sub-surface
    return UnsortedMeshedSurface
    (
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newZones),
        xferCopy(zoneToc_)
    );
}


template<class Face>
Foam::UnsortedMeshedSurface<Face> Foam::UnsortedMeshedSurface<Face>::subsetMesh
(
    const labelHashSet& include
) const
{
    labelList pointMap, faceMap;
    return subsetMesh(include, pointMap, faceMap);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::reset
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst,
    const Xfer<List<label> >& zoneIds
)
{
    ParentType::reset(pointLst, faceLst);

    if (&zoneIds)
    {
        zoneIds_.transfer(zoneIds());
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
        xferMove(surf.zoneIds_)
    );
    zoneToc_.transfer(surf.zoneToc_);

    surf.clear();
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    reset(xferMove(surf.storedPoints()), xferMove(surf.storedFaces()));
    setZones(surf.zones());
    surf.clear();
}


template<class Face>
Foam::Xfer< Foam::UnsortedMeshedSurface<Face> >
Foam::UnsortedMeshedSurface<Face>::xfer()
{
    return xferMove(*this);
}



// Read from file, determine format from extension
template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::read(const fileName& name)
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
bool Foam::UnsortedMeshedSurface<Face>::read
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
    zoneIds_ = surf.zoneIds_;
    zoneToc_ = surf.zoneToc_;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UnsortedMeshedSurfaceIO.C"
#include "UnsortedMeshedSurfaceNew.C"

// ************************************************************************* //

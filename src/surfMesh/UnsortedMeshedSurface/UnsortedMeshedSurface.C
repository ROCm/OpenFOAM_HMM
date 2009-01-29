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
    const Xfer<List<label> >& regionIds,
    const Xfer<surfRegionIdentifierList>& regionTofc
)
:
    ParentType(pointLst, faceLst),
    regionIds_(regionIds),
    regionToc_(regionTofc)
{}


template<class Face>
Foam::UnsortedMeshedSurface<Face>::UnsortedMeshedSurface
(
    const Xfer<pointField>& pointLst,
    const Xfer<List<Face> >& faceLst,
    const UList<label>& regionSizes,
    const UList<word>& regionNames
)
:
    ParentType(pointLst, faceLst)
{
    if (&regionSizes)
    {
        if (&regionNames)
        {
            setRegions(regionSizes, regionNames);
        }
        else
        {
            setRegions(regionSizes);
        }
    }
    else
    {
        oneRegion();
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
    setRegions(surf.regions());
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
    regionIds_(surf.regionIds_),
    regionToc_(surf.regionToc_)
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
void Foam::UnsortedMeshedSurface<Face>::oneRegion(const word& name)
{
    regionIds_.setSize(size());
    regionIds_ = 0;

    word regionName(name);
    if (regionName.empty())
    {
        if (regionToc_.size())
        {
            regionName = regionToc_[0].name();
        }
        if (regionName.empty())
        {
            regionName = "region0";
        }
    }

    // set single default region
    regionToc_.setSize(1);
    regionToc_[0] = surfRegionIdentifier(regionName, 0);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setRegions
(
    const surfRegionList& regionLst
)
{
    regionIds_.setSize(size());
    regionToc_.setSize(regionLst.size());

    forAll(regionToc_, regionI)
    {
        const surfRegion& reg = regionLst[regionI];

        regionToc_[regionI] = reg;

        // assign sub-region Ids
        SubList<label> subRegion(regionIds_, reg.size(), reg.start());
        subRegion = regionI;
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setRegions
(
    const UList<label>& sizes,
    const UList<word>& names
)
{
    regionIds_.setSize(size());
    regionToc_.setSize(sizes.size());

    label start = 0;
    forAll(regionToc_, regionI)
    {
        regionToc_[regionI] = surfRegionIdentifier(names[regionI], regionI);

        // assign sub-region Ids
        SubList<label> subRegion(regionIds_, sizes[regionI], start);
        subRegion = regionI;

        start += sizes[regionI];
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setRegions
(
    const UList<label>& sizes
)
{
    regionIds_.setSize(size());
    regionToc_.setSize(sizes.size());

    label start = 0;
    forAll(regionToc_, regionI)
    {
        regionToc_[regionI] = surfRegionIdentifier
        (
            word("region") + ::Foam::name(regionI),
            regionI
        );

        // assign sub-region Ids
        SubList<label> subRegion(regionIds_, sizes[regionI], start);
        subRegion = regionI;

        start += sizes[regionI];
    }
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::remapFaces
(
    const UList<label>& faceMap
)
{
    // re-assign the region Ids
    if (&faceMap && faceMap.size())
    {
        if (regionToc_.empty())
        {
            oneRegion();
        }
        else if (regionToc_.size() == 1)
        {
            // optimized for single-region case
            regionIds_ = 0;
        }
        else
        {
            List<label> newRegions(faceMap.size());

            forAll(faceMap, faceI)
            {
                newRegions[faceI] = regionIds_[faceMap[faceI]];
            }
            regionIds_.transfer(newRegions);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::UnsortedMeshedSurface<Face>::setSize(const label s)
{
    ParentType::setSize(s);
    // if regions extend: set with last regionId
    regionIds_.setSize(s, regionToc_.size() - 1);
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::clear()
{
    ParentType::clear();
    regionIds_.clear();
    regionToc_.clear();
}


template<class Face>
Foam::surfRegionList Foam::UnsortedMeshedSurface<Face>::sortedRegions
(
    labelList& faceMap
) const
{
    // supply some region names
    Map<word> regionNames;
    forAll(regionToc_, regionI)
    {
        regionNames.insert(regionI, regionToc_[regionI].name());
    }

    return sortedRegionsById(regionIds_, regionNames, faceMap);
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

        newRegions[faceI] = regionIds_[origFaceI];
    }
    oldToNew.clear();

    // construct a sub-surface
    return UnsortedMeshedSurface
    (
        xferMove(newPoints),
        xferMove(newFaces),
        xferMove(newRegions),
        xferCopy(regionToc_)
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
    const Xfer<List<label> >& regionIds
)
{
    ParentType::reset(pointLst, faceLst);

    if (&regionIds)
    {
        regionIds_.transfer(regionIds());
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
        xferMove(surf.regionIds_)
    );
    regionToc_.transfer(surf.regionToc_);

    surf.clear();
}


template<class Face>
void Foam::UnsortedMeshedSurface<Face>::transfer
(
    MeshedSurface<Face>& surf
)
{
    reset(xferMove(surf.storedPoints()), xferMove(surf.storedFaces()));
    setRegions(surf.regions());
    surf.clear();
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
    regionIds_ = surf.regionIds_;
    regionToc_ = surf.regionToc_;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UnsortedMeshedSurfaceIO.C"
#include "UnsortedMeshedSurfaceNew.C"

// ************************************************************************* //

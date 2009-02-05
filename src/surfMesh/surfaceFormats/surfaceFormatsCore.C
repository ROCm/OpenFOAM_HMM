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

#include "surfaceFormatsCore.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
#include "SortableList.H"
#include "surfMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::fileFormats::surfaceFormatsCore::nativeExt("ofs");

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool
Foam::fileFormats::surfaceFormatsCore::isNative(const word& ext)
{
    return (ext == nativeExt);
}


Foam::string
Foam::fileFormats::surfaceFormatsCore::getLineNoComment
(
    IFstream& is
)
{
    string line;
    do
    {
        is.getLine(line);
    }
    while ((line.empty() || line[0] == '#') && is.good());

    return line;
}


Foam::fileName
Foam::fileFormats::surfaceFormatsCore::localMeshFileName(const word& surfName)
{
    const word name(surfName.size() ? surfName : surfaceRegistry::defaultName);

    return fileName
    (
        surfaceRegistry::subInstance
      / name
      / surfMesh::meshSubDir
      / name + "." + nativeExt
    );
}


Foam::fileName
Foam::fileFormats::surfaceFormatsCore::findMeshInstance
(
    const Time& d,
    const word& surfName
)
{
    fileName localName = localMeshFileName(surfName);

    // Search back through the time directories list to find the time
    // closest to and lower than current time

    instantList ts = d.times();
    label instanceI;

    for (instanceI = ts.size()-1; instanceI >= 0; --instanceI)
    {
        if (ts[instanceI].value() <= d.timeOutputValue())
        {
            break;
        }
    }

    // Noting that the current directory has already been searched
    // for mesh data, start searching from the previously stored time directory

    if (instanceI >= 0)
    {
        for (label i = instanceI; i >= 0; --i)
        {
            if (isFile(d.path()/ts[i].name()/localName))
            {
                return ts[i].name();
            }
        }
    }

    return "constant";
}


Foam::fileName
Foam::fileFormats::surfaceFormatsCore::findMeshFile
(
    const Time& d,
    const word& surfName
)
{
    fileName localName = localMeshFileName(surfName);

    // Search back through the time directories list to find the time
    // closest to and lower than current time

    instantList ts = d.times();
    label instanceI;

    for (instanceI = ts.size()-1; instanceI >= 0; --instanceI)
    {
        if (ts[instanceI].value() <= d.timeOutputValue())
        {
            break;
        }
    }

    // Noting that the current directory has already been searched
    // for mesh data, start searching from the previously stored time directory

    if (instanceI >= 0)
    {
        for (label i = instanceI; i >= 0; --i)
        {
            fileName testName(d.path()/ts[i].name()/localName);

            if (isFile(testName))
            {
                return testName;
            }
        }
    }

    // fallback to "constant"
    return d.path()/"constant"/localName;
}


// Returns zone info.
// Sets faceMap to the indexing according to zone numbers.
// Zone numbers start at 0.
Foam::surfZoneList
Foam::fileFormats::surfaceFormatsCore::sortedZonesById
(
    const UList<label>& zoneIds,
    const Map<word>& zoneNames,
    labelList& faceMap
)
{
    // determine sort order according to zone numbers

    // std::sort() really seems to mix up the order.
    // and std::stable_sort() might take too long / too much memory

    // Assuming that we have relatively fewer zones compared to the
    // number of items, just do it ourselves

    // step 1: get zone sizes and store (origId => zoneI)
    Map<label> lookup;
    forAll(zoneIds, faceI)
    {
        const label origId = zoneIds[faceI];

        Map<label>::iterator fnd = lookup.find(origId);
        if (fnd != lookup.end())
        {
            fnd()++;
        }
        else
        {
            lookup.insert(origId, 1);
        }
    }

    // step 2: assign start/size (and name) to the newZones
    // re-use the lookup to map (zoneId => zoneI)
    surfZoneList zoneLst(lookup.size());
    label start = 0;
    label zoneI = 0;
    forAllIter(Map<label>, lookup, iter)
    {
        label origId = iter.key();

        word name;
        Map<word>::const_iterator fnd = zoneNames.find(origId);
        if (fnd != zoneNames.end())
        {
            name = fnd();
        }
        else
        {
            name = word("zone") + ::Foam::name(zoneI);
        }

        zoneLst[zoneI] = surfZone
        (
            name,
            0,           // initialize with zero size
            start,
            zoneI
        );

        // increment the start for the next zone
        // and save the (zoneId => zoneI) mapping
        start += iter();
        iter() = zoneI++;
    }


    // step 3: build the re-ordering
    faceMap.setSize(zoneIds.size());

    forAll(zoneIds, faceI)
    {
        label zoneI = lookup[zoneIds[faceI]];
        faceMap[faceI] =
            zoneLst[zoneI].start() + zoneLst[zoneI].size()++;
    }

    // with reordered faces registered in faceMap
    return zoneLst;
}


bool
Foam::fileFormats::surfaceFormatsCore::checkSupport
(
    const wordHashSet& available,
    const word& ext,
    const bool verbose,
    const word& functionName
)
{
    if (available.found(ext))
    {
        return true;
    }
    else if (verbose)
    {
        wordList toc = available.toc();
        SortableList<word> known(toc.xfer());

        Info<<"Unknown file extension for " << functionName
            << " : " << ext << nl
            <<"Valid types: ( " << nativeExt;
        // compact output:
        forAll(known, i)
        {
            Info<<" " << known[i];
        }
        Info<<" )" << endl;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::surfaceFormatsCore::surfaceFormatsCore()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileFormats::surfaceFormatsCore::~surfaceFormatsCore()
{}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

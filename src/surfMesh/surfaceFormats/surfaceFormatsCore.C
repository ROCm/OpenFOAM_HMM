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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::fileFormats::surfaceFormatsCore::meshSubDir("meshedSurface");
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
Foam::fileFormats::surfaceFormatsCore::findMeshInstance
(
    const Time& d,
    const word& subdirName
)
{
    fileName foamName(d.caseName() + "." + nativeExt);

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
            if (file(d.path()/ts[i].name()/subdirName/foamName))
            {
                return ts[i].name();
            }
        }
    }

    return "constant";
}


Foam::fileName
Foam::fileFormats::surfaceFormatsCore::findMeshName
(
    const Time& d,
    const word& subdirName
)
{
    fileName foamName(d.caseName() + "." + nativeExt);

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
            fileName testName(d.path()/ts[i].name()/subdirName/foamName);

            if (file(testName))
            {
                return testName;
            }
        }
    }

    return d.path()/"constant"/subdirName/foamName;
}


Foam::fileName
Foam::fileFormats::surfaceFormatsCore::findMeshInstance
(
    const Time& d
)
{
    return findMeshInstance(d, meshSubDir);
}


Foam::fileName
Foam::fileFormats::surfaceFormatsCore::findMeshName
(
    const Time& d
)
{
    return findMeshName(d, meshSubDir);
}


// Returns region info.
// Sets faceMap to the indexing according to region numbers.
// Region numbers start at 0.
Foam::surfRegionList
Foam::fileFormats::surfaceFormatsCore::sortedRegionsById
(
    const UList<label>& regionIds,
    const Map<word>& regionNames,
    labelList& faceMap
)
{
    // determine sort order according to region numbers

    // std::sort() really seems to mix up the order.
    // and std::stable_sort() might take too long / too much memory

    // Assuming that we have relatively fewer regions compared to the
    // number of items, just do it ourselves

    // step 1: get region sizes and store (regionId => regionI)
    Map<label> lookup;
    forAll(regionIds, faceI)
    {
        const label regId = regionIds[faceI];

        Map<label>::iterator fnd = lookup.find(regId);
        if (fnd != lookup.end())
        {
            fnd()++;
        }
        else
        {
            lookup.insert(regId, 1);
        }
    }

    // step 2: assign start/size (and name) to the newRegions
    // re-use the lookup to map (regionId => regionI)
    surfRegionList regionLst(lookup.size());
    label start = 0;
    label regionI = 0;
    forAllIter(Map<label>, lookup, iter)
    {
        label regId = iter.key();

        word name;
        Map<word>::const_iterator fnd = regionNames.find(regId);
        if (fnd != regionNames.end())
        {
            name = fnd();
        }
        else
        {
            name = word("region") + ::Foam::name(regionI);
        }

        regionLst[regionI] = surfRegion
        (
            name,
            0,           // initialize with zero size
            start,
            regionI
        );

        // increment the start for the next region
        // and save the (regionId => regionI) mapping
        start += iter();
        iter() = regionI++;
    }


    // step 3: build the re-ordering
    faceMap.setSize(regionIds.size());

    forAll(regionIds, faceI)
    {
        label regionI = lookup[regionIds[faceI]];
        faceMap[faceI] =
            regionLst[regionI].start() + regionLst[regionI].size()++;
    }

    // with reordered faces registered in faceMap
    return regionLst;
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

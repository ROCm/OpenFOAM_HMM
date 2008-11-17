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

#include "surfaceFormatsCore.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::fileFormats::surfaceFormatsCore::meshSubDir("meshedSurface");
Foam::word Foam::fileFormats::surfaceFormatsCore::nativeExt("ofs");

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

//- Check if file extension corresponds to 'native' surface format
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
    while ((line.size() == 0 || line[0] == '#') && is.good());

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
            if (file(d.path()/ts[j].name()/subdirName/foamName))
            {
                return ts[j].name();
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
            fileName testName(d.path()/ts[j].name()/subdirName/foamName);

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

// Returns patch info.
// Sets faceMap to the indexing according to patch numbers.
// Patch numbers start at 0.
Foam::surfGroupList
Foam::fileFormats::surfaceFormatsCore::sortedPatchRegions
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

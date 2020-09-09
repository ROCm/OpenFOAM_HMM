/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "surfaceWriterCaching.H"
#include "ListOps.H"
#include "Fstream.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Compare time values with tolerance
static const equalOp<scalar> equalTimes(ROOTSMALL);

// Use ListOps findLower (with tolerance), to find the location of the next
// time-related index.
// The returned index is always 0 or larger (no negative values).
static label findTimeIndex(const UList<scalar>& list, const scalar val)
{
    label idx =
        findLower
        (
            list,
            val,
            0,
            [](const scalar a, const scalar b)
            {
                return (a < b) && (Foam::mag(b - a) > ROOTSMALL);
            }
        );

    if (idx < 0 || !equalTimes(list[idx], val))
    {
        ++idx;
    }

    return idx;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::writerCaching::writerCaching(const word& cacheFileName)
:
    dictName_(cacheFileName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::surfaceWriters::writerCaching::fieldsDict() const
{
    const dictionary* dictptr = cache_.findDict("fields", keyType::LITERAL);

    if (!dictptr)
    {
        dictptr = &dictionary::null;
    }

    return *dictptr;
}


Foam::dictionary& Foam::surfaceWriters::writerCaching::fieldDict
(
    const word& fieldName
)
{
    return
        cache_
            .subDictOrAdd("fields", keyType::LITERAL)
            .subDictOrAdd(fieldName, keyType::LITERAL);
}


bool Foam::surfaceWriters::writerCaching::remove(const word& fieldName)
{
    dictionary* dictptr = cache_.findDict("fields", keyType::LITERAL);

    if (dictptr)
    {
        return dictptr->remove(fieldName);
    }

    return false;
}


void Foam::surfaceWriters::writerCaching::clear()
{
    times_.clear();
    geoms_.clear();
    cache_.clear();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::surfaceWriters::writerCaching::readPreviousTimes
(
    const fileName& dictFile,
    const scalar timeValue
)
{
    // In 1906 and earlier, the fieldsDict contained "meshes" and "times"
    // entries, each with their own time values.
    // This makes it more difficult to define the exact correspondence
    // between geometry intervals and times.
    //
    // Now track the used geometry intervals as a bitSet.


    // Only called from master
    label timeIndex = 0;
    cache_.clear();

    IFstream is(dictFile);

    if (is.good() && cache_.read(is))
    {
        geoms_.clear();

        cache_.readIfPresent("times", times_);
        timeIndex = findTimeIndex(times_, timeValue);

        labelList geomIndices;
        scalarList meshTimes;

        if (cache_.readIfPresent("geometry", geomIndices))
        {
            // Convert indices to bitSet entries
            geoms_.set(geomIndices);
        }
        else if (cache_.readIfPresent("meshes", meshTimes))
        {
            WarningInFunction
                << nl
                << "Setting geometry timeset information from time values"
                << " (cache from an older OpenFOAM version)." << nl
                << "This may not be fully reliable." << nl
                << nl;

            for (const scalar meshTime : meshTimes)
            {
                const label geomIndex = findTimeIndex(times_, meshTime);
                geoms_.set(geomIndex);
            }
        }

        // Make length consistent with time information.
        // We read/write the indices instead of simply dumping the bitSet.
        // This makes the contents more human readable.
        geoms_.resize(times_.size());
    }

    return timeIndex;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfaceWriters::writerCaching::update
(
    const fileName& baseDir,
    const scalar timeValue,
    const bool geomChanged,
    const word& fieldName,
    const word& fieldType,
    const word& varName
)
{
    const fileName dictFile(baseDir/dictName_);

    bool stateChanged = false;

    const label timeIndex =
    (
        times_.empty()
      ? readPreviousTimes(dictFile, timeValue)
      : findTimeIndex(times_, timeValue)
    );


    // Update stored times list and geometry index

    if (timeIndex < geoms_.size()-1)
    {
        // Clear old content when shrinking
        geoms_.unset(timeIndex);
    }

    // Extend or truncate list
    geoms_.resize(timeIndex+1);
    times_.resize(timeIndex+1, VGREAT);

    if (!equalTimes(times_[timeIndex], timeValue))
    {
        stateChanged = true;
        times_[timeIndex] = timeValue;
    }

    if (geomChanged)
    {
        stateChanged = true;
        geoms_.set(timeIndex);
    }


    // Update time/geometry information in dictionary
    cache_.set("times", times_);
    cache_.set("geometry", geoms_.sortedToc());

    // Debugging, or if needed for older versions:
    //// cache_.set
    //// (
    ////     "meshes",
    ////     IndirectList<scalar>(times_, geoms_.sortedToc())
    //// );

    // Add field information to dictionary
    dictionary& dict = fieldDict(fieldName);

    if (dict.empty())
    {
        stateChanged = true;

        dict.set("type", fieldType);
        if (!varName.empty() && varName != fieldName)
        {
            // Use variable name, if it differs from fieldName
            dict.set("name", varName);
        }
    }

    if (stateChanged)
    {
        OFstream os(dictFile);
        os << "// State file for surface writer output" << nl << nl;
        cache_.write(os, false);

        os << nl << "// End" << nl;
    }

    return stateChanged;
}


// ************************************************************************* //

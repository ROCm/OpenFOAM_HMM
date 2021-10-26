/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "timeSelector.H"
#include "ListOps.H"
#include "argList.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeSelector::timeSelector(const std::string& str)
:
    scalarRanges(str)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::timeSelector::selected(const instant& value) const
{
    return scalarRanges::match(value.value());
}


Foam::List<bool> Foam::timeSelector::selected(const instantList& times) const
{
    List<bool> selectTimes(times.size(), false);

    // Check ranges, avoid false positive on constant/
    forAll(times, timei)
    {
        if (times[timei].name() != "constant" && selected(times[timei]))
        {
            selectTimes[timei] = true;
        }
    }

    // Check specific values
    for (const scalarRange& range : *this)
    {
        if (range.single())
        {
            const scalar target = range.value();

            const label nearestIndex =
                TimePaths::findClosestTimeIndex(times, target);

            // Note could also test if the index is too far away.
            // Eg, for times (0 10 20 30 40) selecting 100 will currently
            // return the closest time (40), but perhaps we should limit that
            // to the last deltaT?

            if (nearestIndex >= 0)
            {
                selectTimes[nearestIndex] = true;
            }
        }
    }

    return selectTimes;
}


Foam::instantList Foam::timeSelector::select(const instantList& times) const
{
    return subset(selected(times), times);
}


void Foam::timeSelector::inplaceSelect(instantList& times) const
{
    inplaceSubset(selected(times), times);
}


void Foam::timeSelector::addOptions
(
    const bool constant,
    const bool withZero
)
{
    if (constant)
    {
        argList::addBoolOption
        (
            "constant",
            "Include the 'constant/' dir in the times list"
        );
    }
    if (withZero)
    {
        argList::addBoolOption
        (
            "withZero",
            "Include the '0/' dir in the times list"
        );
    }
    argList::addBoolOption
    (
        "noZero",
        string("Exclude the '0/' dir from the times list")
      + (
            withZero
          ? ", has precedence over the -withZero option"
          : ""
        )
    );
    argList::addBoolOption
    (
        "latestTime",
        "Select the latest time"
    );
    argList::addOption
    (
        "time",
        "ranges",
        "List of ranges. Eg, ':10,20 40:70 1000:', 'none', etc"
    );
}


Foam::instantList Foam::timeSelector::select
(
    const instantList& times,
    const argList& args,
    const word& constantName
)
{
    if (times.size())
    {
        List<bool> selectTimes(times.size(), true);

        label constantIdx = -1;
        label zeroIdx = -1;
        label latestIdx = -1;

        // Determine locations of constant/ and 0/ directories
        forAll(times, timei)
        {
            if (times[timei].name() == constantName)
            {
                constantIdx = timei;
            }
            else if (times[timei].value() == 0)
            {
                zeroIdx = timei;
            }

            if (constantIdx >= 0 && zeroIdx >= 0)
            {
                break;
            }
        }

        // Determine latestTime selection (if any)
        // This must appear before the -time option processing
        if (args.found("latestTime"))
        {
            selectTimes = false;
            latestIdx = times.size() - 1;

            // Avoid false match on constant/
            if (latestIdx == constantIdx)
            {
                latestIdx = -1;
            }
        }

        if (args.found("time"))
        {
            // Can match 0/, but can never match constant/
            selectTimes = timeSelector(args["time"]).selected(times);
        }

        // Add in latestTime (if selected)
        if (latestIdx >= 0)
        {
            selectTimes[latestIdx] = true;
        }

        if (constantIdx >= 0)
        {
            // Only add constant/ if specifically requested
            selectTimes[constantIdx] = args.found("constant");
        }

        // Special treatment for 0/
        if (zeroIdx >= 0)
        {
            if (args.found("noZero"))
            {
                // Exclude 0/ if specifically requested
                selectTimes[zeroIdx] = false;
            }
            else if (argList::validOptions.found("withZero"))
            {
                // With -withZero enabled, drop 0/ unless specifically requested
                selectTimes[zeroIdx] = args.found("withZero");
            }
        }

        return subset(selectTimes, times);
    }

    return times;
}


Foam::instantList Foam::timeSelector::select0
(
    Time& runTime,
    const argList& args
)
{
    instantList times
    (
        timeSelector::select
        (
            runTime.times(),
            args,
            runTime.constant()
        )
    );

    if (times.empty())
    {
        WarningInFunction
            << "No time specified or available, selecting 'constant'"
            << endl;

        times.append(instant(0, runTime.constant()));
    }

    runTime.setTime(times.first(), 0);

    return times;
}


Foam::instantList Foam::timeSelector::selectIfPresent
(
    Time& runTime,
    const argList& args
)
{
    if
    (
        args.found("latestTime")
     || args.found("time")
     || args.found("constant")
     || args.found("noZero")
     || args.found("withZero")
    )
    {
        return select0(runTime, args);
    }

    // No timeSelector option specified. Do not change runTime.
    return instantList(one{}, instant(runTime.value(), runTime.timeName()));
}


// ************************************************************************* //

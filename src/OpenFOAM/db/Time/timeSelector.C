/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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
#include "StringStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeSelector::timeSelector()
:
    scalarRanges()
{}


Foam::timeSelector::timeSelector(Istream& is)
:
    scalarRanges(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::timeSelector::selected(const instant& value) const
{
    return scalarRanges::selected(value.value());
}


Foam::List<bool> Foam::timeSelector::selected(const instantList& times) const
{
    List<bool> lst(times.size(), false);

    // Check ranges, avoid false positive on constant/
    forAll(times, timei)
    {
        if (times[timei].name() != "constant" && selected(times[timei]))
        {
            lst[timei] = true;
        }
    }

    // Check specific values
    for (const scalarRange& range : *this)
    {
        if (range.isExact())
        {
            const scalar target = range.value();

            int nearestIndex = -1;
            scalar nearestDiff = Foam::GREAT;

            forAll(times, timei)
            {
                if (times[timei].name() == "constant") continue;

                scalar diff = fabs(times[timei].value() - target);
                if (diff < nearestDiff)
                {
                    nearestDiff = diff;
                    nearestIndex = timei;
                }
            }

            if (nearestIndex >= 0)
            {
                lst[nearestIndex] = true;
            }
        }
    }

    return lst;
}


Foam::instantList Foam::timeSelector::select(const instantList& times)
const
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
            "include the 'constant/' dir in the times list"
        );
    }
    if (withZero)
    {
        argList::addBoolOption
        (
            "withZero",
            "include the '0/' dir in the times list"
        );
    }
    argList::addBoolOption
    (
        "noZero",
        string("exclude the '0/' dir from the times list")
      + (
            withZero
          ? ", has precedence over the -withZero option"
          : ""
        )
    );
    argList::addBoolOption
    (
        "latestTime",
        "select the latest time"
    );
    argList::addBoolOption
    (
        "newTimes",
        "select the new times"
    );
    argList::addOption
    (
        "time",
        "ranges",
        "comma-separated time ranges - eg, ':10,20,40:70,1000:'"
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
        if (args.optionFound("latestTime"))
        {
            selectTimes = false;
            latestIdx = times.size() - 1;

            // Avoid false match on constant/
            if (latestIdx == constantIdx)
            {
                latestIdx = -1;
            }
        }

        if (args.optionFound("time"))
        {
            // Can match 0/, but can never match constant/
            selectTimes = timeSelector
            (
                args.optionLookup("time")()
            ).selected(times);
        }

        // Add in latestTime (if selected)
        if (latestIdx >= 0)
        {
            selectTimes[latestIdx] = true;
        }

        if (constantIdx >= 0)
        {
            // Only add constant/ if specifically requested
            selectTimes[constantIdx] = args.optionFound("constant");
        }

        // Special treatment for 0/
        if (zeroIdx >= 0)
        {
            if (args.optionFound("noZero"))
            {
                // Exclude 0/ if specifically requested
                selectTimes[zeroIdx] = false;
            }
            else if (argList::validOptions.found("withZero"))
            {
                // With -withZero enabled, drop 0/ unless specifically requested
                selectTimes[zeroIdx] = args.optionFound("withZero");
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

    runTime.setTime(times[0], 0);

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
        args.optionFound("latestTime")
     || args.optionFound("time")
     || args.optionFound("constant")
     || args.optionFound("noZero")
     || args.optionFound("withZero")
    )
    {
        return select0(runTime, args);
    }

    // No timeSelector option specified. Do not change runTime.
    return instantList(1, instant(runTime.value(), runTime.timeName()));
}


Foam::instantList Foam::timeSelector::select
(
    Time& runTime,
    const argList& args,
    const word& fName
)
{
    instantList times(timeSelector::select0(runTime, args));

    if (times.size() && args.optionFound("newTimes"))
    {
        List<bool> selectTimes(times.size(), true);

        forAll(times, timei)
        {
            selectTimes[timei] =
               !fileHandler().exists
                (
                    runTime.path()
                  / times[timei].name()
                  / fName
                );
        }

        return subset(selectTimes, times);
    }

    return times;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "parProfiling.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "profilingPstream.H"
#include "Tuple2.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(parProfiling, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        parProfiling,
        dictionary
    );

} // End namespace functionObject
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::parProfiling::parProfiling
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name)
{
    profilingPstream::enable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::parProfiling::~parProfiling()
{
    profilingPstream::disable();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::parProfiling::report()
{
    if (!profilingPstream::active())
    {
        return;
    }

    // (Time, Processor) for each of: min/max/sum
    typedef FixedList<Tuple2<double, int>, 3> statData;
    typedef FixedList<statData, 3> statDataTimes;

    // Reduction: if x and y are unequal assign value.
    auto statsEqOp = [](statDataTimes& xStats, const statDataTimes& yStats)
    {
        forAll(xStats, i)
        {
            statData& x = xStats[i];
            const statData& y = yStats[i];

            // 0: min, 1: max, 2: total (or avg)
            if (x[0].first() > y[0].first())
            {
                x[0] = y[0];
            }
            if (x[1].first() < y[1].first())
            {
                x[1] = y[1];
            }
            x[2].first() += y[2].first();
        }
    };

    statDataTimes times;

    // Master time
    {
        const double total =
        (
            profilingPstream::times(profilingPstream::REDUCE)
          + profilingPstream::times(profilingPstream::GATHER)
          + profilingPstream::times(profilingPstream::SCATTER)
            // Include broadcast with reduce instead of all-to-all
          + profilingPstream::times(profilingPstream::BROADCAST)
        );

        times[0] = Tuple2<double, int>(total, Pstream::myProcNo());
    }

    // All time
    {
        const double total =
        (
            profilingPstream::times(profilingPstream::WAIT)
          + profilingPstream::times(profilingPstream::ALL_TO_ALL)
          + profilingPstream::times(profilingPstream::OTHER)
        );

        times[1] = Tuple2<double, int>(total, Pstream::myProcNo());
    }

    // Other time
    {
        const double total =
        (
            profilingPstream::times(profilingPstream::OTHER)
        );

        times[2] = Tuple2<double, int>(total, Pstream::myProcNo());
    }

    profilingPstream::suspend();

    Pstream::combineGather(times, statsEqOp);

    profilingPstream::resume();


    if (Pstream::master())
    {
        Info<< type() << ':' << nl
            << incrIndent;

        {
            const statData& stats = times[0];
            double avg = stats[2].first()/Pstream::nProcs();

            Info<< indent << "reduce    : avg = " << avg << 's' << nl
                << indent << "            min = " << stats[0].first()
                << "s (processor " << stats[0].second() << ')' << nl
                << indent << "            max = " << stats[1].first()
                << "s (processor " << stats[1].second() << ')' << nl;
        }

        {
            const statData& stats = times[1];
            double avg = stats[2].first()/Pstream::nProcs();

            Info<< indent << "all-all   : avg = " << avg << 's' << nl
                << indent << "            min = " << stats[0].first()
                << "s (processor " << stats[0].second() << ')' << nl
                << indent << "            max = " << stats[1].first()
                << "s (processor " << stats[1].second() << ')' << nl;
        }

        {
            const statData& stats = times[2];
            double avg = stats[2].first()/Pstream::nProcs();

            Info<< indent << "other     : avg = " << avg << 's' << nl
                << indent << "            min = " << stats[0].first()
                << "s (processor " << stats[0].second() << ')' << nl
                << indent << "            max = " << stats[1].first()
                << "s (processor " << stats[1].second() << ')' << nl;
        }

        Info<< decrIndent;
    }
}


bool Foam::functionObjects::parProfiling::execute()
{
    report();
    return true;
}


bool Foam::functionObjects::parProfiling::write()
{
    return true;
}


bool Foam::functionObjects::parProfiling::end()
{
    profilingPstream::disable();
    return true;
}


// ************************************************************************* //

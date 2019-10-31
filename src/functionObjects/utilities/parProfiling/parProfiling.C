/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
#include "UPstream.H"
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


    // Processor and time for each of: -min -max -sum
    typedef FixedList<Tuple2<label, scalar>, 3> statData;


    //- Reduction class. If x and y are not equal assign value.
    struct statsEqOp
    {
        void operator()
        (
            FixedList<statData, 2>& xStats,
            const FixedList<statData, 2>& yStats
        ) const
        {
            forAll(xStats, i)
            {
                statData& x = xStats[i];
                const statData& y = yStats[i];

                // 0 : min
                // 1 : max
                // 2 : sum
                if (x[0].second() > y[0].second())
                {
                    x[0].second() = y[0].second();
                    x[0].first()  = y[0].first();
                }
                if (x[1].second() < y[1].second())
                {
                    x[1].second() = y[1].second();
                    x[1].first()  = y[1].first();
                }
                x[2].second() += y[2].second();
                x[2].first()++;
            }
        }
    };

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

    typedef FixedList<Tuple2<label, scalar>, 3> statData;
    FixedList<statData, 2> times;

    {
        const scalar masterTime =
        (
            profilingPstream::times(profilingPstream::REDUCE)
          + profilingPstream::times(profilingPstream::GATHER)
          + profilingPstream::times(profilingPstream::SCATTER)
        );

        statData& reduceStats = times[0];

        Tuple2<label, scalar>& minTime = reduceStats[0];
        minTime.first() = Pstream::myProcNo();
        minTime.second() = masterTime;

        Tuple2<label, scalar>& maxTime = reduceStats[1];
        maxTime.first() = Pstream::myProcNo();
        maxTime.second() = masterTime;

        Tuple2<label, scalar>& sumTime = reduceStats[2];
        sumTime.first() = 1;
        sumTime.second() = masterTime;
    }

    {
        const scalar allTime =
        (
            profilingPstream::times(profilingPstream::WAIT)
          + profilingPstream::times(profilingPstream::ALL_TO_ALL)
        );

        statData& allToAllStats = times[1];

        Tuple2<label, scalar>& minTime = allToAllStats[0];
        minTime.first() = Pstream::myProcNo();
        minTime.second() = allTime;

        Tuple2<label, scalar>& maxTime = allToAllStats[1];
        maxTime.first() = Pstream::myProcNo();
        maxTime.second() = allTime;

        Tuple2<label, scalar>& sumTime = allToAllStats[2];
        sumTime.first() = 1;
        sumTime.second() = allTime;
    }

    profilingPstream::suspend();

    Pstream::combineGather(times, statsEqOp());

    profilingPstream::resume();


    if (Pstream::master())
    {
        const statData& reduceStats = times[0];
        const statData& allToAllStats = times[1];

        scalar reduceAvg = reduceStats[2].second()/Pstream::nProcs();
        scalar allToAllAvg = allToAllStats[2].second()/Pstream::nProcs();

        Info<< type() << ':' << nl
            << incrIndent
            << indent << "reduce    : avg = " << reduceAvg << 's' << nl
            << indent << "            min = " << reduceStats[0].second()
            << "s (processor " << reduceStats[0].first() << ')' << nl
            << indent << "            max = " << reduceStats[1].second()
            << "s (processor " << reduceStats[1].first() << ')' << nl
            << indent << "all-all   : avg = " << allToAllAvg << 's' << nl
            << indent << "            min = " << allToAllStats[0].second()
            << "s (processor " << allToAllStats[0].first() << ')' << nl
            << indent << "            max = " << allToAllStats[1].second()
            << "s (processor " << allToAllStats[1].first() << ')'
            << decrIndent << endl;
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

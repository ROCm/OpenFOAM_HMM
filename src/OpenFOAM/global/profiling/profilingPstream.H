/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::profilingPstream

Description
    Timers and values for simple (simplistic) mpi-profiling. The entire
    class behaves as a singleton.

SourceFiles
    profilingPstream.C

\*---------------------------------------------------------------------------*/

#ifndef profilingPstream_H
#define profilingPstream_H

#include "autoPtr.H"
#include "cpuTime.H"
#include "scalar.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class profilingPstream Declaration
\*---------------------------------------------------------------------------*/

class profilingPstream
{
    //- Timer to use
    static autoPtr<cpuTime> timer_;

    //- Stash for timer during suspend
    static autoPtr<cpuTime> suspend_;

    //- The timing values
    static FixedList<scalar, 5> times_;


public:

    //- Enumeration within times array
    enum timingType
    {
        GATHER = 0,
        SCATTER,
        REDUCE,
        WAIT,
        ALL_TO_ALL
    };

public:

    // Constructors

        //- Construct and enable global timer
        profilingPstream();


    //- Destructor - remove global timer
    ~profilingPstream();


    // Member Functions

        //- Create timer for measuring communication, or reset existing
        static void enable();

        //- Remove timer for measuring communication activity
        static void disable();

        //- Suspend use of timer (if active)
        static void suspend();

        //- Resume use of timer (if previously active)
        static void resume();

        //- Timer is active
        inline static bool active()
        {
            return timer_.valid();
        }

        //- Access to the timing information
        inline static FixedList<scalar, 5>& times()
        {
            return times_;
        }

        //- Access to the timing information
        inline static scalar times(const enum timingType idx)
        {
            return times_[idx];
        }

        //- Update timer prior to measurement
        inline static void beginTiming()
        {
            if (timer_.valid())
            {
                (void) timer_->cpuTimeIncrement();
            }
        }

        //- Add time increment
        inline static void addTime(const enum timingType idx)
        {
            if (timer_.valid())
            {
                times_[idx] += timer_->cpuTimeIncrement();
            }
        }

        //- Add time increment to gatherTime
        inline static void addGatherTime()
        {
            addTime(GATHER);
        }

        //- Add time increment to scatterTime
        inline static void addScatterTime()
        {
            addTime(SCATTER);
        }

        //- Add time increment to reduceTime
        inline static void addReduceTime()
        {
            addTime(REDUCE);
        }

        //- Add time increment to waitTime
        inline static void addWaitTime()
        {
            addTime(WAIT);
        }

        //- Add time increment to allToAllTime
        inline static void addAllToAllTime()
        {
            addTime(ALL_TO_ALL);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

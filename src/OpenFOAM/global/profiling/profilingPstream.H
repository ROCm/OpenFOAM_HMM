/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Timers and values for simple (simplistic) mpi-profiling.
    The entire class behaves as a singleton.

SourceFiles
    profilingPstream.C

\*---------------------------------------------------------------------------*/

#ifndef profilingPstream_H
#define profilingPstream_H

#include "cpuTime.H"
#include "FixedList.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class profilingPstream Declaration
\*---------------------------------------------------------------------------*/

class profilingPstream
{
public:

    // Public Types

        //- Enumeration within times array
        enum timingType
        {
            GATHER = 0,
            SCATTER,
            REDUCE,
            WAIT,
            ALL_TO_ALL
        };

        //- The timing values
        typedef FixedList<double, 5> timingList;


private:

    // Private Static Data

        //- The timer to use
        static std::unique_ptr<cpuTime> timer_;

        //- The timing values
        static timingList times_;

        //- Is timer in a suspend state?
        static bool suspend_;


public:

    // Constructors

        //- Default construct, enables global timer
        profilingPstream();


    //- Destructor, disables global timer
    ~profilingPstream();


    // Member Functions

        //- Create timer for measuring communication, or reset existing
        static void enable();

        //- Remove timer for measuring communication activity
        static void disable();

        //- Suspend use of timer (if active)
        inline static void suspend()
        {
            suspend_ = bool(timer_);
        }

        //- Resume use of timer (if previously active)
        static void resume()
        {
            suspend_ = false;
        }

        //- Timer is active
        inline static bool active()
        {
            return !suspend_ && bool(timer_);
        }

        //- Access to the timing information
        inline static timingList& times()
        {
            return times_;
        }

        //- Access to the timing information at given index
        inline static double times(const enum timingType idx)
        {
            return times_[idx];
        }

        //- Update timer prior to measurement
        inline static void beginTiming()
        {
            if (active())
            {
                (void) timer_->cpuTimeIncrement();
            }
        }

        //- Add time increment
        inline static void addTime(const enum timingType idx)
        {
            if (active())
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

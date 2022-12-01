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

Class
    Foam::profilingPstream

Description
    Timers and values for simple (simplistic) mpi-profiling.
    The entire class behaves as a singleton.

SourceFiles
    profilingPstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_profilingPstream_H
#define Foam_profilingPstream_H

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
            BROADCAST,
            REDUCE,
            WAIT,
            ALL_TO_ALL,
            OTHER
        };

        //- The timing values
        typedef FixedList<double, 7> timingList;


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
        static void disable() noexcept;

        //- Suspend use of timer (if active)
        static void suspend() noexcept
        {
            suspend_ = bool(timer_);
        }

        //- Resume use of timer (if previously active)
        static void resume() noexcept
        {
            suspend_ = false;
        }

        //- Timer is active (not suspended and enabled)
        static bool active() noexcept
        {
            return !suspend_ && bool(timer_);
        }

        //- Access to the timing information
        static timingList& times() noexcept
        {
            return times_;
        }

        //- Access to the timing information at given index
        static double times(const timingType idx)
        {
            return times_[idx];
        }

        //- Update timer prior to measurement
        static void beginTiming()
        {
            if (active())
            {
                (void) timer_->cpuTimeIncrement();
            }
        }

        //- Add time increment
        static void addTime(const timingType idx)
        {
            if (active())
            {
                times_[idx] += timer_->cpuTimeIncrement();
            }
        }

        //- Add time increment to \em gather time
        static void addGatherTime()
        {
            addTime(timingType::GATHER);
        }

        //- Add time increment to \em scatter time
        static void addScatterTime()
        {
            addTime(timingType::SCATTER);
        }

        //- Add time increment to \em broadcast time
        static void addBroadcastTime()
        {
            addTime(timingType::BROADCAST);
        }

        //- Add time increment to \em reduce time
        static void addReduceTime()
        {
            addTime(timingType::REDUCE);
        }

        //- Add time increment to \em wait time
        static void addWaitTime()
        {
            addTime(timingType::WAIT);
        }

        //- Add time increment to \em allToAll time
        static void addAllToAllTime()
        {
            addTime(timingType::ALL_TO_ALL);
        }

        //- Add time increment to \em other time
        static void addOtherTime()
        {
            addTime(timingType::OTHER);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

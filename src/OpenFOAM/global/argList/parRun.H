/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::ParRunControl

Description
    Helper class for initializing parallel jobs from the command arguments,
    storing 'dry-run' state etc.
    Also handles cleanup of parallel (or serial) jobs.

Note
    In the meantime the class name may be slightly misleading.

\*---------------------------------------------------------------------------*/

#ifndef argListRunControl_H
#define argListRunControl_H

#include "Pstream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ParRunControl Declaration
\*---------------------------------------------------------------------------*/

class ParRunControl
{
    int dryRun_;
    int verbose_;
    bool parallel_;
    bool distributed_;

public:

    //- Default construct
    ParRunControl()
    :
        dryRun_(0),
        verbose_(0),
        parallel_(false),
        distributed_(false)
    {}

    //- Destructor. Shutdown (finalize) MPI as required
    ~ParRunControl()
    {
        if (parallel_)
        {
            Info<< "Finalising parallel run" << endl;
        }
        UPstream::shutdown();
    }


    // Parallel Control

        //- Initialize Pstream for a parallel run
        void runPar(int& argc, char**& argv, bool needsThread)
        {
            if (!UPstream::init(argc, argv, needsThread))
            {
                Info<< "Failed to start parallel run" << endl;
                UPstream::exit(1);
            }
            parallel_ = true;
        }


        //- True if this is a parallel run
        bool parRun() const noexcept
        {
            return parallel_;
        }

        //- True if this is a parallel run and uses distributed roots.
        bool distributed() const noexcept
        {
            return (parallel_ && distributed_);
        }

        //- Set use of distributed roots, but only if actually parallel
        void distributed(bool on) noexcept
        {
            distributed_ = (parallel_ && on);
        }


    // General Control

        //- Non-zero if set as 'dry-run'
        int dryRun() const noexcept
        {
            return dryRun_;
        }

        //- Change 'dry-run', return old value
        int dryRun(const int level) noexcept
        {
            int old(dryRun_);
            dryRun_ = level;
            return old;
        }

        //- Non-zero if set as 'verbose'
        int verbose() const noexcept
        {
            return verbose_;
        }

        //- Change 'verbose', return old value
        int verbose(const int level) noexcept
        {
            int old(verbose_);
            verbose_ = level;
            return old;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

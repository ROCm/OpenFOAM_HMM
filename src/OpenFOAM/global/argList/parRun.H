/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

Class
    Foam::ParRunControl

Description
    Helper class for initializing parallel jobs from the command arguments.
    Also handles cleanup of parallel (or serial) jobs.

\*---------------------------------------------------------------------------*/

#ifndef parRun_H
#define parRun_H

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
    bool parallel_;
    bool distributed_;

public:

    //- Default construct
    ParRunControl()
    :
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

        Pstream::shutdown();
    }


    //- Initialize Pstream for a parallel run
    void runPar(int& argc, char**& argv, bool needsThread)
    {
        if (!Pstream::init(argc, argv, needsThread))
        {
            Info<< "Failed to start parallel run" << endl;
            Pstream::exit(1);
        }
        parallel_ = true;
    }

    //- True if this is parallel run.
    bool parRun() const
    {
        return parallel_;
    }

    //- True if this is a parallel run and uses distributed roots.
    bool distributed() const
    {
        return parallel_ && distributed_;
    }

    //- Set use of distributed roots.
    void distributed(bool on)
    {
        distributed_ = (parallel_ ? on : false);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

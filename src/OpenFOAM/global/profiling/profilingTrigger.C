/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009-2016 Bernhard Gschaider
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "profiling.H"
#include "profilingTrigger.H"
#include "profilingInformation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profilingTrigger::profilingTrigger()
:
    ptr_(nullptr)
{}


Foam::profilingTrigger::profilingTrigger(const char* name)
:
    ptr_(profiling::New(name))
{}


Foam::profilingTrigger::profilingTrigger(const string& name)
:
    ptr_(profiling::New(name))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profilingTrigger::~profilingTrigger()
{
    stop();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::profilingTrigger::running() const
{
    return ptr_;
}


void Foam::profilingTrigger::stop()
{
    if (ptr_)
    {
        // profiling info pointer managed by pool storage, so no delete here
        profiling::unstack(ptr_);
    }

    ptr_ = nullptr;
}


// ************************************************************************* //

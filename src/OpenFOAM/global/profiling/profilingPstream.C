/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "profilingPstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

std::unique_ptr<Foam::cpuTime> Foam::profilingPstream::timer_(nullptr);

bool Foam::profilingPstream::suspend_(false);

Foam::profilingPstream::timingList Foam::profilingPstream::times_(double(0));
Foam::profilingPstream::countList Foam::profilingPstream::counts_(uint64_t(0));


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profilingPstream::profilingPstream()
{
    enable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profilingPstream::~profilingPstream()
{
    disable();
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::profilingPstream::enable()
{
    if (timer_)
    {
        timer_->resetCpuTime(); // Not necessarily required ...
    }
    else
    {
        timer_.reset(new cpuTime);
        times_ = double(0);
        counts_ = uint64_t(0);
    }

    suspend_ = false;
}


void Foam::profilingPstream::disable() noexcept
{
    timer_.reset(nullptr);
    suspend_ = false;
}


double Foam::profilingPstream::elapsedTime()
{
    double total = 0;
    for (const double val : times_)
    {
        total += val;
    }

    return total;
}


// ************************************************************************* //

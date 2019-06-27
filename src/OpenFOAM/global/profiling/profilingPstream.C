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

\*---------------------------------------------------------------------------*/

#include "profilingPstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cpuTime> Foam::profilingPstream::timer_(nullptr);

Foam::autoPtr<Foam::cpuTime> Foam::profilingPstream::suspend_(nullptr);

Foam::FixedList<Foam::scalar, 5> Foam::profilingPstream::times_(Zero);


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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::profilingPstream::enable()
{
    if (timer_.valid())
    {
        timer_->resetCpuTime(); // Not really needed ...
    }
    else if (suspend_.valid())
    {
        suspend_.swap(timer_);
        timer_->resetCpuTime(); // Not really needed ...
    }
    else
    {
        timer_.reset(new cpuTime);
        times_ = Zero;
    }

    suspend_.clear();
}


void Foam::profilingPstream::disable()
{
    timer_.clear();
    suspend_.clear();
}


void Foam::profilingPstream::suspend()
{
    suspend_.clear();
    suspend_.swap(timer_);
}


void Foam::profilingPstream::resume()
{
    if (suspend_.valid())
    {
        timer_.clear();
        timer_.swap(suspend_);
    }
}


// ************************************************************************* //

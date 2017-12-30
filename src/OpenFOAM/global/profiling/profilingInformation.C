/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2016 Bernhard Gschaider
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "profilingInformation.H"
#include "Switch.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::profilingInformation::nextId_(0);


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::profilingInformation::getNextId()
{
    return nextId_++;
}


void Foam::profilingInformation::raiseId(label maxVal)
{
    if (nextId_ < maxVal)
    {
        nextId_ = maxVal;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profilingInformation::profilingInformation()
:
    id_(getNextId()),
    description_("application::main"),
    parent_(this),
    calls_(0),
    totalTime_(0),
    childTime_(0),
    maxMem_(0),
    onStack_(false)
{}


Foam::profilingInformation::profilingInformation
(
    const string& descr,
    profilingInformation *parent
)
:
    id_(getNextId()),
    description_(descr),
    parent_(parent),
    calls_(0),
    totalTime_(0),
    childTime_(0),
    maxMem_(0),
    onStack_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::profilingInformation::update(const scalar elapsed)
{
    ++calls_;
    totalTime_ += elapsed;

    if (id_ != parent().id())
    {
        parent().childTime_ += elapsed;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::profilingInformation::push() const
{
    onStack_ = true;
}


void Foam::profilingInformation::pop() const
{
    onStack_ = false;
}


Foam::Ostream& Foam::profilingInformation::write
(
    Ostream& os,
    const bool offset,
    const scalar elapsedTime,
    const scalar childTimes
) const
{
    // write in dictionary format

    os.beginBlock(word("trigger" + Foam::name(id_)));

    os.writeEntry("id",             id_);
    os.writeEntryIfDifferent("parentId", id_, parent().id());
    os.writeEntry("description",    description());
    os.writeEntry("calls",          calls()     + (offset ? 1 : 0));
    os.writeEntry("totalTime",      totalTime() + elapsedTime);
    os.writeEntry("childTime",      childTime() + childTimes);
    os.writeEntryIfDifferent<int>("maxMem", 0, maxMem_);
    os.writeEntry("onStack",        Switch(onStack()));

    os.endBlock();

    return os;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const profilingInformation& info)
{
    return info.write(os);
}


// ************************************************************************* //

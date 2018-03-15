/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2016 Bernhard Gschaider
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "argList.H"
#include "profiling.H"
#include "profilingInformation.H"
#include "profilingSysInfo.H"
#include "cpuInfo.H"
#include "memInfo.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::profiling::allowed
(
    Foam::debug::infoSwitch("allowProfiling", 1)
);

Foam::profiling* Foam::profiling::pool_(nullptr);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::profilingInformation* Foam::profiling::find
(
    const string& descr,
    const label parentId
)
{
    return hash_.lookup(Key(descr, parentId), nullptr);
}


Foam::profilingInformation* Foam::profiling::store(profilingInformation *info)
{
    // Profile information lookup is qualified by parent id
    hash_.insert(Key(info->description(), info->parent().id()), info);
    return info;
}


void Foam::profiling::push(profilingInformation *info, clockTime& timer)
{
    stack_.push(info);
    timers_.set(info->id(), &timer);
    info->push();                       // mark as on stack
}


Foam::profilingInformation* Foam::profiling::pop()
{
    profilingInformation *info = stack_.pop();
    timers_.erase(info->id());
    info->pop();                        // mark as off stack

    return info;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::profiling::active()
{
    return allowed && pool_;
}


void Foam::profiling::disable()
{
    allowed = 0;
}


bool Foam::profiling::print(Ostream& os)
{
    if (active())
    {
        return pool_->writeData(os);
    }

    return false;
}


bool Foam::profiling::writeNow()
{
    if (active())
    {
        return pool_->regIOobject::write();
    }

    return false;
}


void Foam::profiling::initialize
(
    const IOobject& ioObj,
    const Time& owner
)
{
    if (allowed && !pool_)
    {
        pool_ = new profiling(ioObj, owner);

        profilingInformation *info = pool_->store
        (
            new profilingInformation()
        );

        pool_->push(info, pool_->clockTime_);
        if (argList::bannerEnabled())
        {
            Info<< "profiling initialized" << nl;
        }
    }

    // silently ignore multiple initializations
    // eg, decomposePar uses multiple simultaneous Times
}


void Foam::profiling::initialize
(
    const dictionary& dict,
    const IOobject& ioObj,
    const Time& owner
)
{
    if (allowed && !pool_)
    {
        pool_ = new profiling(dict, ioObj, owner);

        profilingInformation *info = pool_->store
        (
            new profilingInformation()
        );

        pool_->push(info, pool_->clockTime_);
        if (argList::bannerEnabled())
        {
            Info<< "profiling initialized" << nl;
        }
    }

    // silently ignore multiple initializations
    // eg, decomposePar uses multiple simultaneous Times
}


void Foam::profiling::stop(const Time& owner)
{
    if (pool_ && &owner == &(pool_->owner_))
    {
        delete pool_;
        pool_ = nullptr;
    }
}


Foam::profilingInformation* Foam::profiling::New
(
    const string& descr,
    clockTime& timer
)
{
    profilingInformation *info = nullptr;

    if (active())
    {
        profilingInformation *parent = pool_->stack_.top();

        info = pool_->find(descr, parent->id());
        if (!info)
        {
            info = pool_->store(new profilingInformation(descr, parent));
        }

        pool_->push(info, timer);

        if (pool_->memInfo_)
        {
            info->maxMem_ = Foam::max
            (
                info->maxMem_,
                pool_->memInfo_->update().size()
            );
        }
    }

    return info;
}


void Foam::profiling::unstack(const profilingInformation *info)
{
    if (active() && info)
    {
        profilingInformation *top = pool_->pop();

        if (info->id() != top->id())
        {
            FatalErrorInFunction
                << "The profiling information to unstack has different"
                << " id than on the top of the profiling stack" << nl
                << "  info: " << info->id() << " (" << info->description()
                << ")\n"
                << "  top:  " << top->id()  << " (" << top->description()
                << ")\n" << endl
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profiling::profiling
(
    const IOobject& io,
    const Time& owner
)
:
    IOdictionary(io),
    owner_(owner),
    clockTime_(),
    hash_(),
    stack_(),
    timers_(),
    sysInfo_(new profilingSysInfo()),
    cpuInfo_(new cpuInfo()),
    memInfo_(new memInfo())
{}


Foam::profiling::profiling
(
    const dictionary& dict,
    const IOobject& io,
    const Time& owner
)
:
    IOdictionary(io),
    owner_(owner),
    clockTime_(),
    hash_(),
    stack_(),
    timers_(),
    sysInfo_
    (
        dict.lookupOrDefault<bool>("sysInfo", false)
      ? new profilingSysInfo() : nullptr
    ),
    cpuInfo_
    (
        dict.lookupOrDefault<bool>("cpuInfo", false)
      ? new cpuInfo() : nullptr
    ),
    memInfo_
    (
        dict.lookupOrDefault<bool>("memInfo", false)
      ? new memInfo() : nullptr
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profiling::~profiling()
{
    deleteDemandDrivenData(sysInfo_);
    deleteDemandDrivenData(cpuInfo_);
    deleteDemandDrivenData(memInfo_);

    if (pool_ == this)
    {
        pool_ = nullptr;
        profilingInformation::nextId_ = 0;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Time& Foam::profiling::owner() const
{
    return owner_;
}


Foam::label Foam::profiling::size() const
{
    return stack_.size();
}


bool Foam::profiling::writeData(Ostream& os) const
{
    os.beginBlock("profiling");

    // Add extra new line between entries
    label nTrigger = 0;

    // write on-stack items
    // newest is first on the stack, top-level is at the end
    // this is how the child times are summed
    {
        scalar oldElapsed = 0;
        forAllConstIter(StackContainer, stack_, iter)
        {
            const profilingInformation *info = *iter;
            scalar elapsed = timers_[info->id()]->elapsedTime();

            if (nTrigger++)
            {
                os << nl;
            }
            info->write(os, true, elapsed, oldElapsed);
            oldElapsed = elapsed;
        }
    }


    // write off-stack items
    // using an additional Map to sort by Id
    {
        typedef Map<const Information*> LookupContainer;
        LookupContainer lookup;

        forAllConstIter(StorageContainer, hash_, iter)
        {
            const profilingInformation *info = iter();

            if (!info->onStack())
            {
                lookup.set(info->id(), info);
            }
        }

        forAllConstIter(LookupContainer, lookup, iter)
        {
            if (nTrigger++)
            {
                os << nl;
            }
            iter()->write(os);
        }
    }

    os.endBlock();

    if (sysInfo_)
    {
        os << nl;
        os.beginBlock("sysInfo");
        sysInfo_->write(os);
        os.endBlock();
    }

    if (cpuInfo_)
    {
        os << nl;
        os.beginBlock("cpuInfo");
        cpuInfo_->write(os);
        os.endBlock();
    }

    if (memInfo_)
    {
        memInfo_->update();

        os << nl;
        os.beginBlock("memInfo");
        memInfo_->write(os);
        os.writeEntry("units", "kB");
        os.endBlock();
    }

    return os;
}


bool Foam::profiling::writeObject
(
    IOstream::streamFormat,
    IOstream::versionNumber ver,
    IOstream::compressionType,
    const bool valid
) const
{
    return regIOobject::writeObject
    (
        IOstream::ASCII,
        ver,
        IOstream::UNCOMPRESSED,
        true
    );
}


// ************************************************************************* //

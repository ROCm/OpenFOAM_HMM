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

#include "Profiling.H"
#include "ProfilingSysInfo.H"
#include "cpuInfo.H"
#include "memInfo.H"
#include "OSspecific.H"
#include "IOstreams.H"
#include "dictionary.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::Profiling* Foam::Profiling::pool_(0);

Foam::label Foam::Profiling::Information::nextId_(0);


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// file-scope function
template<class T>
inline static void writeEntry
(
    Foam::Ostream& os, const Foam::word& key, const T& value
)
{
    os.writeKeyword(key) << value << Foam::token::END_STATEMENT << '\n';
}


Foam::label Foam::Profiling::Information::getNextId()
{
    return nextId_++;
}


void Foam::Profiling::Information::raiseID(label maxVal)
{
    if (nextId_ < maxVal)
    {
        nextId_ = maxVal;
    }
}


bool Foam::Profiling::active()
{
    return pool_;
}


bool Foam::Profiling::writeNow()
{
    if (pool_)
    {
        Info<<"Profiling::writeNow() at time = "
            << pool_->owner().timeName() << endl;
        return pool_->write();
    }
    else
    {
        return false;
    }
}


void Foam::Profiling::initialize
(
    const IOobject& ioObj,
    const Time& owner
)
{
    if (pool_)
    {
        WarningInFunction
            << "Already initialized" << endl;
    }
    else
    {
        pool_ = new Profiling(ioObj, owner);

        Information *info = pool_->store
        (
            new Information()
        );

        pool_->push(info, pool_->clockTime_);
        Info<< "Profiling initialized" << nl;
    }
}


void Foam::Profiling::initialize
(
    const dictionary& dict,
    const IOobject& ioObj,
    const Time& owner
)
{
    if (pool_)
    {
        WarningInFunction
            << "Already initialized" << endl;
    }
    else
    {
        pool_ = new Profiling(dict, ioObj, owner);

        Information *info = pool_->store
        (
            new Information()
        );

        pool_->push(info, pool_->clockTime_);
        Info<< "Profiling initialized" << nl;
    }
}


void Foam::Profiling::stop(const Time& owner)
{
    if (pool_ && &owner == &(pool_->owner_))
    {
        delete pool_;
        pool_ = 0;
    }
}


Foam::Profiling::Information* Foam::Profiling::New
(
    const string& name,
    clockTime& timer
)
{
    Information *info = 0;

    if (pool_)
    {
        info = pool_->find(name);
        if (!info)
        {
            info = pool_->store
            (
                new Information(pool_->stack_.top(), name)
            );
        }

        pool_->push(info, timer);
    }

    return info;
}


void Foam::Profiling::unstack(const Information *info)
{
    if (pool_ && info)
    {
        Information *top = pool_->pop();

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

Foam::Profiling::Profiling
(
    const IOobject& io,
    const Time& owner
)
:
    regIOobject(io),
    owner_(owner),
    clockTime_(),
    hash_(),
    stack_(),
    timers_(),
    sysInfo_(new sysInfo()),
    cpuInfo_(new cpuInfo()),
    memInfo_(new memInfo())
{}


Foam::Profiling::Profiling
(
    const dictionary& dict,
    const IOobject& io,
    const Time& owner
)
:
    regIOobject(io),
    owner_(owner),
    clockTime_(),
    hash_(),
    stack_(),
    timers_(),
    sysInfo_
    (
        dict.lookupOrDefault<Switch>("sysInfo", true)
      ? new sysInfo() : 0
    ),
    cpuInfo_
    (
        dict.lookupOrDefault<Switch>("cpuInfo", true)
      ? new cpuInfo() : 0
    ),
    memInfo_
    (
        dict.lookupOrDefault<Switch>("memInfo", true)
      ? new memInfo() : 0
    )
{}


Foam::Profiling::Information::Information()
:
    id_(getNextId()),
    description_("application::main"),
    parent_(this),
    calls_(0),
    totalTime_(0),
    childTime_(0),
    onStack_(false)
{}


Foam::Profiling::Information::Information
(
    Information *parent,
    const string& descr
)
:
    id_(getNextId()),
    description_(descr),
    parent_(parent),
    calls_(0),
    totalTime_(0),
    childTime_(0),
    onStack_(false)
{}


Foam::Profiling::Trigger::Trigger(const char* name)
:
    clock_(),
    ptr_(Profiling::New(name, clock_))
{}


Foam::Profiling::Trigger::Trigger(const string& name)
:
    clock_(),
    ptr_(Profiling::New(name, clock_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Profiling::~Profiling()
{
    deleteDemandDrivenData(sysInfo_);
    deleteDemandDrivenData(cpuInfo_);
    deleteDemandDrivenData(memInfo_);

    if (pool_ == this)
    {
        pool_ = 0;
        Information::nextId_ = 0;
    }
}


Foam::Profiling::Information::~Information()
{}


Foam::Profiling::Trigger::~Trigger()
{
    stop();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Time& Foam::Profiling::owner() const
{
    return owner_;
}

Foam::label Foam::Profiling::size() const
{
    return stack_.size();
}


Foam::Profiling::Information* Foam::Profiling::find(const string& name)
{
    StorageContainer::iterator iter = hash_.find(name);
    return (iter != hash_.end() ? iter() : 0);
}


void Foam::Profiling::Information::update(const scalar& elapsed)
{
    ++calls_;
    totalTime_ += elapsed;

    if (id_ != parent().id())
    {
        parent().childTime_ += elapsed;
    }
}


bool Foam::Profiling::writeData(Ostream& os) const
{
    os  << indent << "profiling" << nl
        << indent << token::BEGIN_LIST << incrIndent << nl;

    // write on-stack items
    // newest is first on the stack, top-level is at the end
    // this is how the child times are summed
    {
        scalar oldElapsed = 0;
        forAllConstIter(StackContainer, stack_, iter)
        {
            const Information *info = *iter;
            scalar elapsed = timers_[info->id()]->elapsedTime();

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
            const Information *info = iter();

            if (!info->onStack())
            {
                lookup.set(info->id(), info);
            }
        }

        forAllConstIter(LookupContainer, lookup, iter)
        {
            iter()->write(os);
        }
    }

    os  << decrIndent
        << indent << token::END_LIST << token::END_STATEMENT << nl;


    if (sysInfo_)
    {
        os << nl;
        os.beginBlock("sysInfo") << nl; // FUTURE: without nl
        sysInfo_->write(os);
        os.endBlock() << nl; // FUTURE: without nl
    }

    if (cpuInfo_)
    {
        os << nl;
        os.beginBlock("cpuInfo") << nl; // FUTURE: without nl
        cpuInfo_->write(os);
        os.endBlock() << nl; // FUTURE: without nl
    }

    if (memInfo_)
    {
        memInfo_->update();

        os << nl;
        os.beginBlock("memInfo") << nl; // FUTURE: without nl
        memInfo_->write(os);
        writeEntry(os, "units", "kB");
        os.endBlock() << nl; // FUTURE: without nl
    }

    return os;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Profiling::Information* Foam::Profiling::store(Information *info)
{
    hash_.insert(info->description(), info);
    return info;
}


void Foam::Profiling::push(Information *info, clockTime& timer)
{
    stack_.push(info);
    timers_.set(info->id(), &timer);
    info->push();                       // mark as on stack
}


Foam::Profiling::Information* Foam::Profiling::pop()
{
    Information *info = stack_.pop();
    timers_.erase(info->id());
    info->pop();                        // mark as off stack

    return info;
}


bool Foam::Profiling::Trigger::running() const
{
    return ptr_;
}


void Foam::Profiling::Trigger::stop()
{
    if (ptr_)
    {
        ptr_->update(clock_.elapsedTime());
        Profiling::unstack(ptr_);
        // pointer is managed by pool storage -> thus no delete here
    }
    ptr_ = 0;
}


void Foam::Profiling::Information::push() const
{
    onStack_ = true;
}


void Foam::Profiling::Information::pop() const
{
    onStack_ = false;
}


Foam::Ostream& Foam::Profiling::Information::write
(
    Ostream& os,
    const bool offset,
    const scalar& elapsedTime,
    const scalar& childTimes
) const
{
    // write in dictionary format

    // os.beginBlock("_" + Foam::name(id_)) << nl;
    os.beginBlock() << nl; // FUTURE: without nl

    // FUTURE: os.writeEntry(key, value);

    writeEntry(os, "id",            id_);
    if (id_ != parent().id())
    {
        writeEntry(os, "parentId",  parent().id());
    }
    writeEntry(os, "description",   description());
    writeEntry(os, "calls",         calls()     + (offset ? 1 : 0));
    writeEntry(os, "totalTime",     totalTime() + elapsedTime);
    writeEntry(os, "childTime",     childTime() + childTimes);
    writeEntry(os, "onStack",       Switch(onStack()));

    os.endBlock() << nl; // FUTURE: without nl

    return os;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Profiling::Information& info
)
{
    return info.write(os);
}


// ************************************************************************* //

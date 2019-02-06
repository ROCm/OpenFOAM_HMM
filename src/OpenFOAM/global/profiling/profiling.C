/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2009-2016 Bernhard Gschaider
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

int Foam::profiling::allowed(Foam::debug::infoSwitch("allowProfiling", 1));
Foam::profiling* Foam::profiling::singleton_(nullptr);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::profilingInformation* Foam::profiling::create(const zero)
{
    // Top-level entry: reset everything
    pool_.clear();
    children_.clear();
    stack_.clear();
    times_.clear();

    Information* info = new Information();

    pool_.append(info);
    children_.resize(pool_.size());
    children_.last().clear();  // safety

    return info;
}


Foam::profilingInformation* Foam::profiling::create
(
    profilingInformation *parent,
    const string& descr
)
{
    const label parentId = parent->id();

    for (Information* child : children_[parentId])
    {
        if (descr == child->description())
        {
            return child;  // Found existing
        }
    }

    Information* info = new Information(parent, descr, pool_.size());

    pool_.append(info);
    children_.resize(pool_.size());
    children_.last().clear();  // safety
    children_[parentId].append(info);

    return info;
}


void Foam::profiling::beginTimer(profilingInformation *info)
{
    stack_.append(info);
    times_.append(clockValue::now());
    info->setActive(true);              // Mark as on stack
}


Foam::profilingInformation* Foam::profiling::endTimer()
{
    Information *info = stack_.remove();
    clockValue clockval = times_.remove();

    info->update(clockval.elapsed());   // Update elapsed time
    info->setActive(false);             // Mark as off stack

    return info;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::profiling::active()
{
    return allowed && singleton_;
}


void Foam::profiling::disable()
{
    allowed = 0;
}


bool Foam::profiling::print(Ostream& os)
{
    if (active())
    {
        return singleton_->writeData(os);
    }

    return false;
}


bool Foam::profiling::writeNow()
{
    if (active())
    {
        return singleton_->regIOobject::write();
    }

    return false;
}


void Foam::profiling::initialize
(
    const IOobject& ioObj,
    const Time& owner
)
{
    if (allowed && !singleton_)
    {
        singleton_ = new profiling(ioObj, owner);
    }
}


void Foam::profiling::initialize
(
    const dictionary& dict,
    const IOobject& ioObj,
    const Time& owner
)
{
    if (allowed && !singleton_)
    {
        singleton_ = new profiling(dict, ioObj, owner);
    }
}


void Foam::profiling::stop(const Time& owner)
{
    if (singleton_ && &owner == &(singleton_->owner_))
    {
        delete singleton_;
        singleton_ = nullptr;
    }
}


Foam::profilingInformation* Foam::profiling::New(const string& descr)
{
    Information *info = nullptr;

    if (active())
    {
        Information *parent = singleton_->stack_.last();

        info = singleton_->create(parent, descr);
        singleton_->beginTimer(info);

        if (singleton_->memInfo_)
        {
            info->maxMem_ = Foam::max
            (
                info->maxMem_,
                singleton_->memInfo_->update().size()
            );
        }
    }

    return info;
}


void Foam::profiling::unstack(const profilingInformation *info)
{
    if (active() && info)
    {
        Information *top = singleton_->endTimer();

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
    pool_(),
    children_(),
    stack_(),
    times_(),
    sysInfo_(new profilingSysInfo()),
    cpuInfo_(new cpuInfo()),
    memInfo_(new memInfo())
{
    Information *info = this->create(Zero);
    this->beginTimer(info);

    DetailInfo << "profiling initialized" << nl;
}


Foam::profiling::profiling
(
    const dictionary& dict,
    const IOobject& io,
    const Time& owner
)
:
    IOdictionary(io),
    owner_(owner),
    pool_(),
    children_(),
    stack_(),
    times_(),
    sysInfo_
    (
        dict.lookupOrDefault("sysInfo", false)
      ? new profilingSysInfo() : nullptr
    ),
    cpuInfo_
    (
        dict.lookupOrDefault("cpuInfo", false)
      ? new cpuInfo() : nullptr
    ),
    memInfo_
    (
        dict.lookupOrDefault("memInfo", false)
      ? new memInfo() : nullptr
    )
{
    Information *info = this->create(Zero);
    this->beginTimer(info);

    DetailInfo << "profiling initialized" << nl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::profiling::~profiling()
{
    deleteDemandDrivenData(sysInfo_);
    deleteDemandDrivenData(cpuInfo_);
    deleteDemandDrivenData(memInfo_);

    if (singleton_ == this)
    {
        singleton_ = nullptr;
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
    static DynamicList<scalar> elapsed;

    const clockValue now(clockValue::now());

    const label nstack = stack_.size();

    elapsed.resize(nstack+1);   // extend for last entry, which has no child.

    for (label stacki=0; stacki < nstack; ++stacki)
    {
        elapsed[stacki] = (now - times_[stacki]);
    }
    elapsed.last() = 0;

    os.beginBlock("profiling");

    // Active items
    for (label stacki=0; stacki < nstack; ++stacki)
    {
        if (stacki) os << nl;   // Extra line between entries

        stack_[stacki]->write
        (
            os,
            true,
            elapsed[stacki],    // elapsedTime
            elapsed[stacki+1]   // childTimes
        );
    }

    // Non-active items
    for (const Information& info : pool_)
    {
        if (!info.active())
        {
            os << nl;
            info.write(os);
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

    return os.good();
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

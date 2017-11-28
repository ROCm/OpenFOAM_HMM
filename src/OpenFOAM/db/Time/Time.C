/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "Time.H"
#include "PstreamReduceOps.H"
#include "argList.H"
#include "HashSet.H"
#include "profiling.H"
#include "demandDrivenData.H"
#include "IOdictionary.H"

#include <sstream>

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Time, 0);
}

const Foam::Enum
<
    Foam::Time::stopAtControls
>
Foam::Time::stopAtControlNames
{
    { stopAtControls::saEndTime, "endTime" },
    { stopAtControls::saNoWriteNow, "noWriteNow" },
    { stopAtControls::saWriteNow, "writeNow" },
    { stopAtControls::saNextWrite, "nextWrite" },
    // NOTE: stopAtControls::saUnknown is left untabulated here so that it can
    // be used as fallback value to flag unknown settings
};


const Foam::Enum
<
    Foam::Time::writeControls
>
Foam::Time::writeControlNames
{
    { writeControls::wcTimeStep, "timeStep" },
    { writeControls::wcRunTime, "runTime" },
    { writeControls::wcAdjustableRunTime, "adjustableRunTime" },
    { writeControls::wcClockTime, "clockTime" },
    { writeControls::wcCpuTime, "cpuTime" },
};


Foam::Time::fmtflags Foam::Time::format_(Foam::Time::general);

int Foam::Time::precision_(6);

const int Foam::Time::maxPrecision_(3 - log10(SMALL));

Foam::word Foam::Time::controlDictName("controlDict");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::Time::adjustDeltaT()
{
    bool adjustTime = false;
    scalar timeToNextWrite = VGREAT;

    if (writeControl_ == wcAdjustableRunTime)
    {
        adjustTime = true;
        timeToNextWrite = max
        (
            0.0,
            (writeTimeIndex_ + 1)*writeInterval_ - (value() - startTime_)
        );
    }

    if (adjustTime)
    {
        scalar nSteps = timeToNextWrite/deltaT_;

        // For tiny deltaT the label can overflow!
        if (nSteps < labelMax)
        {
            // nSteps can be < 1 so make sure at least 1
            label nStepsToNextWrite = max(1, round(nSteps));

            scalar newDeltaT = timeToNextWrite/nStepsToNextWrite;

            // Control the increase of the time step to within a factor of 2
            // and the decrease within a factor of 5.
            if (newDeltaT >= deltaT_)
            {
                deltaT_ = min(newDeltaT, 2.0*deltaT_);
            }
            else
            {
                deltaT_ = max(newDeltaT, 0.2*deltaT_);
            }
        }
    }

    functionObjects_.adjustTimeStep();
}


void Foam::Time::setControls()
{
    // default is to resume calculation from "latestTime"
    const word startFrom = controlDict_.lookupOrDefault<word>
    (
        "startFrom",
        "latestTime"
    );

    if (startFrom == "startTime")
    {
        controlDict_.lookup("startTime") >> startTime_;
    }
    else
    {
        // Search directory for valid time directories
        instantList timeDirs = findTimes(path(), constant());

        if (startFrom == "firstTime")
        {
            if (timeDirs.size())
            {
                if (timeDirs[0].name() == constant() && timeDirs.size() >= 2)
                {
                    startTime_ = timeDirs[1].value();
                }
                else
                {
                    startTime_ = timeDirs[0].value();
                }
            }
        }
        else if (startFrom == "latestTime")
        {
            if (timeDirs.size())
            {
                startTime_ = timeDirs.last().value();
            }
        }
        else
        {
            FatalIOErrorInFunction(controlDict_)
                << "expected startTime, firstTime or latestTime"
                << " found '" << startFrom << "'"
                << exit(FatalIOError);
        }
    }

    setTime(startTime_, 0);

    readDict();
    deltaTSave_ = deltaT_;
    deltaT0_ = deltaT_;

    // Check if time directory exists
    // If not increase time precision to see if it is formatted differently.
    if (!fileHandler().exists(timePath(), false))
    {
        int oldPrecision = precision_;
        int requiredPrecision = -1;
        bool found = false;
        word oldTime(timeName());
        for
        (
            precision_ = maxPrecision_;
            precision_ > oldPrecision;
            precision_--
        )
        {
            // Update the time formatting
            setTime(startTime_, 0);

            word newTime(timeName());
            if (newTime == oldTime)
            {
                break;
            }
            oldTime = newTime;

            // Check the existence of the time directory with the new format
            found = fileHandler().exists(timePath(), false);

            if (found)
            {
                requiredPrecision = precision_;
            }
        }

        if (requiredPrecision > 0)
        {
            // Update the time precision
            precision_ = requiredPrecision;

            // Update the time formatting
            setTime(startTime_, 0);

            WarningInFunction
                << "Increasing the timePrecision from " << oldPrecision
                << " to " << precision_
                << " to support the formatting of the current time directory "
                << timeName() << nl << endl;
        }
        else
        {
            // Could not find time directory so assume it is not present
            precision_ = oldPrecision;

            // Revert the time formatting
            setTime(startTime_, 0);
        }
    }

    if (Pstream::parRun())
    {
        scalar sumStartTime = startTime_;
        reduce(sumStartTime, sumOp<scalar>());
        if
        (
            mag(Pstream::nProcs()*startTime_ - sumStartTime)
          > Pstream::nProcs()*deltaT_/10.0
        )
        {
            FatalIOErrorInFunction(controlDict_)
                << "Start time is not the same for all processors" << nl
                << "processor " << Pstream::myProcNo() << " has startTime "
                << startTime_ << exit(FatalIOError);
        }
    }

    IOdictionary timeDict
    (
        IOobject
        (
            "time",
            timeName(),
            "uniform",
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    // Read and set the deltaT only if time-step adjustment is active
    // otherwise use the deltaT from the controlDict
    if (controlDict_.lookupOrDefault<Switch>("adjustTimeStep", false))
    {
        if (timeDict.readIfPresent("deltaT", deltaT_))
        {
            deltaTSave_ = deltaT_;
            deltaT0_ = deltaT_;
        }
    }

    timeDict.readIfPresent("deltaT0", deltaT0_);

    if (timeDict.readIfPresent("index", startTimeIndex_))
    {
        timeIndex_ = startTimeIndex_;
    }


    // Check if values stored in time dictionary are consistent

    // 1. Based on time name
    bool checkValue = true;

    string storedTimeName;
    if (timeDict.readIfPresent("name", storedTimeName))
    {
        if (storedTimeName == timeName())
        {
            // Same time. No need to check stored value
            checkValue = false;
        }
    }

    // 2. Based on time value
    //    (consistent up to the current time writing precision so it won't
    //     trigger if we just change the write precision)
    if (checkValue)
    {
        scalar storedTimeValue;
        if (timeDict.readIfPresent("value", storedTimeValue))
        {
            word storedTimeName(timeName(storedTimeValue));

            if (storedTimeName != timeName())
            {
                IOWarningInFunction(timeDict)
                    << "Time read from time dictionary " << storedTimeName
                    << " differs from actual time " << timeName() << '.' << nl
                    << "    This may cause unexpected database behaviour."
                    << " If you are not interested" << nl
                    << "    in preserving time state delete"
                    << " the time dictionary."
                    << endl;
            }
        }
    }
}


void Foam::Time::setMonitoring(const bool forceProfiling)
{
    const dictionary* profilingDict = controlDict_.subDictPtr("profiling");
    if (!profilingDict)
    {
        // ... or from etc/controlDict
        profilingDict = debug::controlDict().subDictPtr("profiling");
    }

    // initialize profiling on request
    // otherwise rely on profiling entry within controlDict
    // and skip if 'active' keyword is explicitly set to false
    if (forceProfiling)
    {
        profiling::initialize
        (
            IOobject
            (
                "profiling",
                timeName(),
                "uniform",
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );
    }
    else if
    (
        profilingDict
     && profilingDict->lookupOrDefault<bool>("active", true)
    )
    {
        profiling::initialize
        (
            *profilingDict,
            IOobject
            (
                "profiling",
                timeName(),
                "uniform",
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );
    }

    // Time objects not registered so do like objectRegistry::checkIn ourselves.
    if (runTimeModifiable_)
    {
        // Monitor all files that controlDict depends on
        fileHandler().addWatches(controlDict_, controlDict_.files());
    }

    // Clear dependent files - not needed now
    controlDict_.files().clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Time::Time
(
    const word& controlDictName,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

    loopProfiling_(nullptr),
    libs_(),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(GREAT),
    purgeWrite_(0),
    subCycling_(0),
    writeOnce_(false),
    sigWriteNow_(true, *this),
    sigStopAtWriteNow_(true, *this),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),

    functionObjects_(*this, enableFunctionObjects)
{
    libs_.open(controlDict_, "libs");

    // Explicitly set read flags on objectRegistry so anything constructed
    // from it reads as well (e.g. fvSolution).
    readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    setControls();
    setMonitoring();
}


Foam::Time::Time
(
    const word& controlDictName,
    const argList& args,
    const word& systemName,
    const word& constantName
)
:
    TimePaths
    (
        args.parRunControl().parRun(),
        args.rootPath(),
        args.distributed(),
        args.globalCaseName(),
        args.caseName(),
        systemName,
        constantName
    ),

    objectRegistry(*this),

    loopProfiling_(nullptr),
    libs_(),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(GREAT),
    purgeWrite_(0),
    subCycling_(0),
    writeOnce_(false),
    sigWriteNow_(true, *this),
    sigStopAtWriteNow_(true, *this),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),

    functionObjects_
    (
        *this,
        argList::validOptions.found("withFunctionObjects")
      ? args.optionFound("withFunctionObjects")
      : argList::validOptions.found("noFunctionObjects")
      ? !args.optionFound("noFunctionObjects")
      : false
    )
{
    libs_.open(controlDict_, "libs");

    // Explicitly set read flags on objectRegistry so anything constructed
    // from it reads as well (e.g. fvSolution).
    readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    setControls();

    // '-profiling' = force profiling, ignore controlDict entry
    setMonitoring(args.optionFound("profiling"));
}


Foam::Time::Time
(
    const dictionary& dict,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

    loopProfiling_(nullptr),
    libs_(),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        dict
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(GREAT),
    purgeWrite_(0),
    subCycling_(0),
    writeOnce_(false),
    sigWriteNow_(true, *this),
    sigStopAtWriteNow_(true, *this),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),

    functionObjects_(*this, enableFunctionObjects)
{
    libs_.open(controlDict_, "libs");

    // Explicitly set read flags on objectRegistry so anything constructed
    // from it reads as well (e.g. fvSolution).
    readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    // Since could not construct regIOobject with setting:
    controlDict_.readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    setControls();
    setMonitoring();
}


Foam::Time::Time
(
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

    loopProfiling_(nullptr),
    libs_(),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(GREAT),
    purgeWrite_(0),
    subCycling_(0),
    writeOnce_(false),
    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),

    functionObjects_(*this, enableFunctionObjects)
{
    libs_.open(controlDict_, "libs");
    setMonitoring(); // for profiling etc
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Time::~Time()
{
    deleteDemandDrivenData(loopProfiling_);

    forAllReverse(controlDict_.watchIndices(), i)
    {
        fileHandler().removeWatch(controlDict_.watchIndices()[i]);
    }

    // Destroy function objects first
    functionObjects_.clear();

    // Clean up profiling
    profiling::stop(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::Time::timeName(const scalar t, const int precision)
{
    std::ostringstream buf;
    buf.setf(ios_base::fmtflags(format_), ios_base::floatfield);
    buf.precision(precision);
    buf << t;
    return buf.str();
}


Foam::word Foam::Time::timeName() const
{
    return dimensionedScalar::name();
}


Foam::instantList Foam::Time::times() const
{
    return findTimes(path(), constant());
}


Foam::word Foam::Time::findInstancePath
(
    const fileName& directory,
    const instant& t
) const
{
    // Simplified version: use findTimes (readDir + sort). The expensive
    // bit is the readDir, not the sorting. Tbd: avoid calling findInstancePath
    // from filePath.

    instantList timeDirs = findTimes(path(), constant());
    // Note:
    // - times will include constant (with value 0) as first element.
    //   For backwards compatibility make sure to find 0 in preference
    //   to constant.
    // - list is sorted so could use binary search

    forAllReverse(timeDirs, i)
    {
        if (t.equal(timeDirs[i].value()))
        {
            return timeDirs[i].name();
        }
    }

    return word::null;
}


Foam::word Foam::Time::findInstancePath(const instant& t) const
{
    return findInstancePath(path(), t);
}


Foam::instant Foam::Time::findClosestTime(const scalar t) const
{
    instantList timeDirs = findTimes(path(), constant());

    // There is only one time (likely "constant") so return it
    if (timeDirs.size() == 1)
    {
        return timeDirs[0];
    }

    if (t < timeDirs[1].value())
    {
        return timeDirs[1];
    }
    else if (t > timeDirs.last().value())
    {
        return timeDirs.last();
    }

    label nearestIndex = -1;
    scalar deltaT = GREAT;

    for (label timei=1; timei < timeDirs.size(); ++timei)
    {
        scalar diff = mag(timeDirs[timei].value() - t);
        if (diff < deltaT)
        {
            deltaT = diff;
            nearestIndex = timei;
        }
    }

    return timeDirs[nearestIndex];
}


Foam::label Foam::Time::findClosestTimeIndex
(
    const instantList& timeDirs,
    const scalar t,
    const word& constantName
)
{
    label nearestIndex = -1;
    scalar deltaT = GREAT;

    forAll(timeDirs, timei)
    {
        if (timeDirs[timei].name() == constantName) continue;

        scalar diff = mag(timeDirs[timei].value() - t);
        if (diff < deltaT)
        {
            deltaT = diff;
            nearestIndex = timei;
        }
    }

    return nearestIndex;
}


Foam::label Foam::Time::startTimeIndex() const
{
    return startTimeIndex_;
}


Foam::dimensionedScalar Foam::Time::startTime() const
{
    return dimensionedScalar("startTime", dimTime, startTime_);
}


Foam::dimensionedScalar Foam::Time::endTime() const
{
    return dimensionedScalar("endTime", dimTime, endTime_);
}


Foam::Time::stopAtControls Foam::Time::stopAt() const
{
    return stopAt_;
}


bool Foam::Time::run() const
{
    deleteDemandDrivenData(loopProfiling_);

    bool isRunning = value() < (endTime_ - 0.5*deltaT_);

    if (!subCycling_)
    {
        // Only execute when the condition is no longer true
        // ie, when exiting the control loop
        if (!isRunning && timeIndex_ != startTimeIndex_)
        {
            // Ensure functionObjects execute on last time step
            // (and hence write uptodate functionObjectProperties)
            addProfiling(foExec, "functionObjects.execute()");
            functionObjects_.execute();
            endProfiling(foExec);

            addProfiling(foEnd, "functionObjects.end()");
            functionObjects_.end();
            endProfiling(foEnd);
        }
    }

    if (isRunning)
    {
        if (!subCycling_)
        {
            const_cast<Time&>(*this).readModifiedObjects();

            if (timeIndex_ == startTimeIndex_)
            {
                addProfiling(functionObjects, "functionObjects.start()");
                functionObjects_.start();
            }
            else
            {
                addProfiling(functionObjects, "functionObjects.execute()");
                functionObjects_.execute();
            }

            // Check if the execution of functionObjects require re-reading
            // any files. This moves effect of e.g. 'timeActivatedFileUpdate'
            // one time step forward. Note that we cannot call
            // readModifiedObjects from within timeActivatedFileUpdate since
            // it might re-read the functionObjects themselves (and delete
            // the timeActivatedFileUpdate one)
            if (functionObjects_.filesModified())
            {
                const_cast<Time&>(*this).readModifiedObjects();
            }
        }

        // Update the "is-running" status following the
        // possible side-effects from functionObjects
        isRunning = value() < (endTime_ - 0.5*deltaT_);

        // (re)trigger profiling
        if (profiling::active())
        {
            loopProfiling_ =
                new profilingTrigger("time.run() " + objectRegistry::name());
        }
    }

    return isRunning;
}


bool Foam::Time::loop()
{
    const bool isRunning = run();

    if (isRunning)
    {
        operator++();
    }

    return isRunning;
}


bool Foam::Time::end() const
{
    return value() > (endTime_ + 0.5*deltaT_);
}


bool Foam::Time::stopAt(const stopAtControls stopCtrl) const
{
    if (stopCtrl == stopAtControls::saUnknown)
    {
        return false;
    }

    const bool changed = (stopAt_ != stopCtrl);
    stopAt_ = stopCtrl;
    endTime_ = GREAT;

    // Adjust endTime
    if (stopCtrl == stopAtControls::saEndTime)
    {
        controlDict_.lookup("endTime") >> endTime_;
    }

    return changed;
}


void Foam::Time::setTime(const Time& t)
{
    value() = t.value();
    dimensionedScalar::name() = t.dimensionedScalar::name();
    timeIndex_ = t.timeIndex_;
    fileHandler().setTime(*this);
}


void Foam::Time::setTime(const instant& inst, const label newIndex)
{
    value() = inst.value();
    dimensionedScalar::name() = inst.name();
    timeIndex_ = newIndex;

    IOdictionary timeDict
    (
        IOobject
        (
            "time",
            timeName(),
            "uniform",
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    timeDict.readIfPresent("deltaT", deltaT_);
    timeDict.readIfPresent("deltaT0", deltaT0_);
    timeDict.readIfPresent("index", timeIndex_);
    fileHandler().setTime(*this);
}


void Foam::Time::setTime(const dimensionedScalar& newTime, const label newIndex)
{
    setTime(newTime.value(), newIndex);
}


void Foam::Time::setTime(const scalar newTime, const label newIndex)
{
    value() = newTime;
    dimensionedScalar::name() = timeName(timeToUserTime(newTime));
    timeIndex_ = newIndex;
    fileHandler().setTime(*this);
}


void Foam::Time::setEndTime(const dimensionedScalar& endTime)
{
    setEndTime(endTime.value());
}


void Foam::Time::setEndTime(const scalar endTime)
{
    endTime_ = endTime;
}


void Foam::Time::setDeltaT
(
    const dimensionedScalar& deltaT,
    const bool adjust
)
{
    setDeltaT(deltaT.value(), adjust);
}


void Foam::Time::setDeltaT(const scalar deltaT, const bool adjust)
{
    deltaT_ = deltaT;
    deltaTchanged_ = true;

    if (adjust)
    {
        adjustDeltaT();
    }
}


Foam::TimeState Foam::Time::subCycle(const label nSubCycles)
{
    prevTimeState_.set(new TimeState(*this));  // Fatal if already set

    setTime(*this - deltaT(), (timeIndex() - 1)*nSubCycles);
    deltaT_ /= nSubCycles;
    deltaT0_ /= nSubCycles;
    deltaTSave_ = deltaT0_;

    subCycling_ = nSubCycles;

    return prevTimeState();
}


void Foam::Time::subCycleIndex(const label index)
{
    // Only permit adjustment if sub-cycling was already active
    // and if the index is valid (positive, non-zero).
    // This avoids potential mixups for deleting.

    if (subCycling_ && index > 0)
    {
        subCycling_ = index;
    }
}


void Foam::Time::endSubCycle()
{
    if (subCycling_)
    {
        TimeState::operator=(prevTimeState());
        prevTimeState_.clear();
    }

    subCycling_ = 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::Time& Foam::Time::operator+=(const dimensionedScalar& deltaT)
{
    return operator+=(deltaT.value());
}


Foam::Time& Foam::Time::operator+=(const scalar deltaT)
{
    setDeltaT(deltaT);
    return operator++();
}


Foam::Time& Foam::Time::operator++()
{
    deltaT0_ = deltaTSave_;
    deltaTSave_ = deltaT_;

    // Save old time value and name
    const scalar oldTimeValue = timeToUserTime(value());
    const word oldTimeName = dimensionedScalar::name();

    // Increment time
    setTime(value() + deltaT_, timeIndex_ + 1);

    if (!subCycling_)
    {
        // If the time is very close to zero reset to zero
        if (mag(value()) < 10*SMALL*deltaT_)
        {
            setTime(0.0, timeIndex_);
        }

        if (sigStopAtWriteNow_.active() || sigWriteNow_.active())
        {
            // A signal might have been sent on one processor only
            // Reduce so all decide the same.

            label flag = 0;
            if (sigStopAtWriteNow_.active() && stopAt_ == saWriteNow)
            {
                flag += 1;
            }
            if (sigWriteNow_.active() && writeOnce_)
            {
                flag += 2;
            }
            reduce(flag, maxOp<label>());

            if (flag & 1)
            {
                stopAt_ = saWriteNow;
            }
            if (flag & 2)
            {
                writeOnce_ = true;
            }
        }

        writeTime_ = false;

        switch (writeControl_)
        {
            case wcTimeStep:
                writeTime_ = !(timeIndex_ % label(writeInterval_));
            break;

            case wcRunTime:
            case wcAdjustableRunTime:
            {
                const label writeIndex = label
                (
                    ((value() - startTime_) + 0.5*deltaT_)
                  / writeInterval_
                );

                if (writeIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = writeIndex;
                }
            }
            break;

            case wcCpuTime:
            {
                const label writeIndex = label
                (
                    returnReduce(elapsedCpuTime(), maxOp<double>())
                  / writeInterval_
                );
                if (writeIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = writeIndex;
                }
            }
            break;

            case wcClockTime:
            {
                const label writeIndex = label
                (
                    returnReduce(label(elapsedClockTime()), maxOp<label>())
                  / writeInterval_
                );
                if (writeIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = writeIndex;
                }
            }
            break;
        }


        // Check if endTime needs adjustment to stop at the next run()/end()
        if (!end())
        {
            if (stopAt_ == saNoWriteNow)
            {
                endTime_ = value();
            }
            else if (stopAt_ == saWriteNow)
            {
                endTime_ = value();
                writeTime_ = true;
            }
            else if (stopAt_ == saNextWrite && writeTime_ == true)
            {
                endTime_ = value();
            }
        }

        // Override writeTime if one-shot writing
        if (writeOnce_)
        {
            writeTime_ = true;
            writeOnce_ = false;
        }

        // Adjust the precision of the time directory name if necessary
        if (writeTime_)
        {
            // Tolerance used when testing time equivalence
            const scalar timeTol =
                max(min(pow(10.0, -precision_), 0.1*deltaT_), SMALL);

            // User-time equivalent of deltaT
            const scalar userDeltaT = timeToUserTime(deltaT_);

            // Time value obtained by reading timeName
            scalar timeNameValue = -VGREAT;

            // Check that new time representation differs from old one
            // reinterpretation of the word
            if
            (
                readScalar(dimensionedScalar::name(), timeNameValue)
             && (mag(timeNameValue - oldTimeValue - userDeltaT) > timeTol)
            )
            {
                int oldPrecision = precision_;
                while
                (
                    precision_ < maxPrecision_
                 && readScalar(dimensionedScalar::name(), timeNameValue)
                 && (mag(timeNameValue - oldTimeValue - userDeltaT) > timeTol)
                )
                {
                    precision_++;
                    setTime(value(), timeIndex());
                }

                if (precision_ != oldPrecision)
                {
                    WarningInFunction
                        << "Increased the timePrecision from " << oldPrecision
                        << " to " << precision_
                        << " to distinguish between timeNames at time "
                        << dimensionedScalar::name()
                        << endl;

                    if (precision_ == maxPrecision_)
                    {
                        // Reached maxPrecision limit
                        WarningInFunction
                            << "Current time name " << dimensionedScalar::name()
                            << nl
                            << "    The maximum time precision has been reached"
                               " which might result in overwriting previous"
                               " results."
                            << endl;
                    }

                    // Check if round-off error caused time-reversal
                    scalar oldTimeNameValue = -VGREAT;
                    if
                    (
                        readScalar(oldTimeName, oldTimeNameValue)
                     && (
                            sign(timeNameValue - oldTimeNameValue)
                         != sign(deltaT_)
                        )
                    )
                    {
                        WarningInFunction
                            << "Current time name " << dimensionedScalar::name()
                            << " is set to an instance prior to the "
                               "previous one "
                            << oldTimeName << nl
                            << "    This might result in temporal "
                               "discontinuities."
                            << endl;
                    }
                }
            }
        }
    }

    return *this;
}


Foam::Time& Foam::Time::operator++(int)
{
    return operator++();
}


// ************************************************************************* //

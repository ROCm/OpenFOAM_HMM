/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "loopControl.H"
#include "fvSolution.H"
#include "wordRes.H"
#include "solutionControl.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::loopControl::clear()
{
    total_ = 0;
    interval_ = 0;

    convergenceDict_.clear();
    onLoop_.clear();
    onConverged_.clear();
    onEnd_.clear();

    converged_ = false;
}


void Foam::loopControl::read(const dictionary& dict)
{
    clear();

    bool enabled = dict.lookupOrDefault("enabled", true);

    if (enabled)
    {
        scalar timeStart;
        if (dict.readIfPresent("timeStart", timeStart))
        {
            timeStart = time_.userTimeToTime(timeStart);

            enabled =
            (
                enabled
             && time_.value() >= (timeStart - 0.5*time_.deltaTValue())
            );
        }

        scalar timeEnd;
        if (dict.readIfPresent("timeEnd", timeEnd))
        {
            timeEnd = time_.userTimeToTime(timeEnd);

            enabled =
            (
                enabled
             && time_.value() <= (timeEnd + 0.5*time_.deltaTValue())
            );
        }
    }

    if (!enabled)
    {
        return;
    }

    dict.readIfPresent("iterations", total_);
    dict.readIfPresent("interval", interval_);

    convergenceDict_ = dict.subOrEmptyDict("convergence");

    dict.readIfPresent("onLoop", onLoop_);
    dict.readIfPresent("onConverged", onConverged_);
    dict.readIfPresent("onEnd", onEnd_);
}


bool Foam::loopControl::checkConverged() const
{
    if (convergenceDict_.empty())
    {
        return false;
    }

    HashTable<const fvMesh*> meshes = time_.lookupClass<const fvMesh>();

    bool achieved = true;
    bool checked = false; // safety that some checks were indeed performed

    forAllConstIters(meshes, meshIter)
    {
        const fvMesh& regionMesh = *(meshIter.object());

        const dictionary& solverDict = regionMesh.solverPerformanceDict();

        forAllConstIters(solverDict, iter)
        {
            const entry& dataDictEntry = *iter;

            const word& variableName = dataDictEntry.keyword();

            const scalar absTol =
                convergenceDict_.lookupOrDefault<scalar>(variableName, -1);

            if (absTol > 0)
            {
                // Treat like a SIMPLE control

                Pair<scalar> residuals =
                    solutionControl::maxResidual
                    (
                        regionMesh,
                        dataDictEntry
                    );

                checked = true;
                achieved = achieved && (residuals.first() < absTol);
            }
        }
    }

    return checked && achieved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::loopControl::loopControl
(
    Time& runTime,
    const label nCycles,
    const word& loopName
)
:
    subLoopTime(runTime, nCycles),
    name_(loopName),
    interval_(0),
    convergenceDict_(),
    onLoop_(),
    onConverged_(),
    onEnd_(),
    converged_(false)
{}


Foam::loopControl::loopControl
(
    Time& runTime,
    const dictionary& algorithmDict,
    const word& dictName
)
:
    loopControl(runTime, 0, dictName)
{
    // The loop sub-dictionary
    const dictionary* dictptr = algorithmDict.subDictPtr(dictName);

    if (dictptr)
    {
        // Info<< dictName << *dictptr << endl;
        read(*dictptr);
    }
}


Foam::loopControl::loopControl
(
    Time& runTime,
    const word& algorithmName,
    const word& dictName
)
:
    loopControl(runTime, 0, dictName)
{
    fvSolution fvsol(time_);

    // Eg, PIMPLE or SIMPLE from <system/fvSolution>
    const dictionary* dictptr =
        fvsol.solutionDict().subDictPtr(algorithmName);

    if (dictptr)
    {
        // The loop sub-dictionary
        dictptr = dictptr->subDictPtr(dictName);

        if (dictptr)
        {
            // Info<< dictName << *dictptr << endl;
            read(*dictptr);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::loopControl::~loopControl()
{
    stop();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::loopControl::loop()
{
    bool active = (index_ < total_);   // as per status()

    if (active)
    {
        operator++();

        converged_ = checkConverged();

        if (converged_)
        {
            time_.functionObjects().execute(onConverged_, index_);
            stop();
            return false;
        }
        else if
        (
            interval_ && !(index_ % interval_)
         && !onLoop_.empty()
        )
        {
            time_.functionObjects().execute(onLoop_, index_);
        }
    }
    else if (index_)
    {
        // Not active, the loop condition has now exiting on the last subloop

        if (!converged_ && !onEnd_.empty())
        {
            time_.functionObjects().execute(onEnd_, index_);
        }
    }

    return active;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const loopControl& ctrl)
{
    os << ctrl.name() << ": ";
    if (ctrl.nCycles() && ctrl.index() <= ctrl.nCycles())
    {
        os << ctrl.index() << '/' << ctrl.nCycles();
    }
    else
    {
        os << "off";
    }

    return os;
}


// ************************************************************************* //

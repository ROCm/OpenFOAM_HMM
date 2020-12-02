/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "SIBaseIncompressible.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(SIBase, 0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void SIBase::read()
{
    surfaceSensitivity_.read();
    includeObjective_ =
        dict().getOrDefault<bool>("includeObjectiveContribution", true);
    writeSensitivityMap_ =
        dict().getOrDefault<bool>("writeSensitivityMap", false);

    // If includeObjective is set to true both here and in the surface
    // sensitivities, set the one in the latter to false to avoid double
    // contributions
    bool surfSensIncludeObjective(surfaceSensitivity_.getIncludeObjective());
    if (includeObjective_ && surfSensIncludeObjective)
    {
        WarningInFunction
            << "includeObjectiveContribution set to true in both "
            << "surfaceSensitivities and the parameterization options" << nl
            << "This will lead to double contributions " << nl
            << "Disabling the former"
            << endl;
        surfaceSensitivity_.setIncludeObjective(false);
    }

    // Make sure surface area is included in the sensitivity map
    surfaceSensitivity_.setIncludeSurfaceArea(true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SIBase::SIBase
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleVars& primalVars,
    incompressibleAdjointVars& adjointVars,
    objectiveManager& objectiveManager
)
:
    shapeSensitivities
    (
        mesh,
        dict,
        primalVars,
        adjointVars,
        objectiveManager
    ),
    surfaceSensitivity_
    (
        mesh,
        // Ideally, subOrEmptyDict would be used.
        // Since we need a recursive search in shapeSensitivities though
        // and the dict returned by subOrEmptyDict (if found)
        // does not know its parent, optionalSubDict is used
        dict.optionalSubDict("surfaceSensitivities"),
        primalVars,
        adjointVars,
        objectiveManager
    ),
    includeObjective_(true),
    writeSensitivityMap_(true)
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool SIBase::readDict(const dictionary& dict)
{
    if (sensitivity::readDict(dict))
    {
        surfaceSensitivity_.readDict
        (
            dict.optionalSubDict("surfaceSensitivities")
        );

        return true;
    }

    return false;
}


void SIBase::accumulateIntegrand(const scalar dt)
{
    // Accumulate multiplier of dxFace/db
    surfaceSensitivity_.accumulateIntegrand(dt);

    // Accumulate direct sensitivities
    if (includeObjective_)
    {
        accumulateDirectSensitivityIntegrand(dt);
    }

    // Accumulate sensitivities due to boundary conditions
    accumulateBCSensitivityIntegrand(dt);
}


void SIBase::clearSensitivities()
{
    surfaceSensitivity_.clearSensitivities();
    shapeSensitivities::clearSensitivities();
}


const sensitivitySurface& SIBase::getSurfaceSensitivities() const
{
    return surfaceSensitivity_;
}


void SIBase::write(const word& baseName)
{
    shapeSensitivities::write(baseName);
    if (writeSensitivityMap_)
    {
        surfaceSensitivity_.write(baseName);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

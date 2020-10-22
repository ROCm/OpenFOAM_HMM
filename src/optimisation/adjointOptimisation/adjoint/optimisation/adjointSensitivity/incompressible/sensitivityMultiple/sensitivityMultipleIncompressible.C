/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "sensitivityMultipleIncompressible.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivityMultiple, 0);
addToRunTimeSelectionTable
(
    adjointSensitivity,
    sensitivityMultiple,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivityMultiple::sensitivityMultiple
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleVars& primalVars,
    incompressibleAdjointVars& adjointVars,
    objectiveManager& objectiveManager
)
:
    adjointSensitivity
    (
        mesh,
        dict,
        primalVars,
        adjointVars,
        objectiveManager
    ),
    sensTypes_(dict.subDict("sensTypes").toc()),
    sens_(sensTypes_.size())
{
    forAll(sensTypes_, sI)
    {
        sens_.set
        (
            sI,
            adjointSensitivity::New
            (
                mesh,
                dict.subDict("sensTypes").subDict(sensTypes_[sI]),
                primalVars,
                adjointVars,
                objectiveManager
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool sensitivityMultiple::readDict(const dictionary& dict)
{
    if (sensitivity::readDict(dict))
    {
        forAll(sens_, sI)
        {
            sens_[sI].readDict
            (
                dict.subDict("sensTypes").subDict(sensTypes_[sI])
            );
        }

        return true;
    }

    return false;
}


void sensitivityMultiple::accumulateIntegrand(const scalar dt)
{
    forAll(sens_, sI)
    {
        sens_[sI].accumulateIntegrand(dt);
    }
}


void sensitivityMultiple::assembleSensitivities()
{
    forAll(sens_, sI)
    {
        sens_[sI].assembleSensitivities();
    }
}


const scalarField& sensitivityMultiple::calculateSensitivities()
{
    forAll(sens_, sI)
    {
        Info<< "Computing sensitivities " << sensTypes_[sI] << endl;
        derivatives_ = sens_[sI].calculateSensitivities();
    }
    write(type());

    return derivatives_;
}


void sensitivityMultiple::clearSensitivities()
{
    forAll(sens_, sI)
    {
        sens_[sI].clearSensitivities();
    }
}


void sensitivityMultiple::write(const word& baseName)
{
    forAll(sens_, sI)
    {
        sens_[sI].write(sensTypes_[sI]);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "incompressibleAdjointVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(incompressibleAdjointVars, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

incompressibleAdjointVars::incompressibleAdjointVars
(
    fvMesh& mesh,
    solverControl& SolverControl,
    objectiveManager& objManager,
    incompressibleVars& primalVars
)
:
    incompressibleAdjointMeanFlowVars(mesh, SolverControl, primalVars),
    objectiveManager_(objManager),

    adjointTurbulence_
    (
        incompressibleAdjoint::adjointRASModel::New
        (
            primalVars_,
            *this,
            objManager
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const autoPtr<incompressibleAdjoint::adjointRASModel>&
incompressibleAdjointVars::adjointTurbulence() const
{
    return adjointTurbulence_;
}


autoPtr<incompressibleAdjoint::adjointRASModel>&
incompressibleAdjointVars::adjointTurbulence()
{
    return adjointTurbulence_;
}


void incompressibleAdjointVars::resetMeanFields()
{
    if (solverControl_.average())
    {
        Info<< "Resetting adjoint mean fields to zero" << endl;

        // Reset fields to zero
        paMeanPtr_() == dimensionedScalar(paPtr_().dimensions(), Zero);
        UaMeanPtr_() == dimensionedVector(UaPtr_().dimensions(), Zero);
        phiaMeanPtr_() == dimensionedScalar(phiaPtr_().dimensions(), Zero);
        adjointTurbulence_().resetMeanFields();

        // Reset averaging iteration index to 0
        solverControl_.averageIter() = 0;
    }
}


void incompressibleAdjointVars::computeMeanFields()
{
    if (solverControl_.doAverageIter())
    {
        Info<< "Averaging adjoint fields" << endl;
        label& iAverageIter = solverControl_.averageIter();
        scalar avIter(iAverageIter);
        scalar oneOverItP1 = 1./(avIter+1);
        scalar mult = avIter*oneOverItP1;
        paMeanPtr_() == paMeanPtr_()  *mult + paPtr_()  *oneOverItP1;
        UaMeanPtr_() == UaMeanPtr_()  *mult + UaPtr_()  *oneOverItP1;
        phiaMeanPtr_() == phiaMeanPtr_()*mult + phiaPtr_()*oneOverItP1;
        adjointTurbulence_().computeMeanFields();
        ++iAverageIter;
    }
}


void incompressibleAdjointVars::nullify()
{
    incompressibleAdjointMeanFlowVars::nullify();
    adjointTurbulence_->nullify();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

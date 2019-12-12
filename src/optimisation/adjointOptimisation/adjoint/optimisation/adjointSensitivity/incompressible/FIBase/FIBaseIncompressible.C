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

#include "FIBaseIncompressible.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(FIBase, 0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void FIBase::read()
{
    includeDistance_ =
        dict_.lookupOrDefault<bool>
        (
            "includeDistance",
            adjointVars_.adjointTurbulence().ref().includeDistance()
        );

    // Allocate distance solver if needed
    if (includeDistance_ && eikonalSolver_.empty())
    {
        eikonalSolver_.reset
        (
            new adjointEikonalSolver
            (
                mesh_,
                dict_,
                primalVars_.RASModelVariables(),
                adjointVars_.adjointTurbulence(),
                sensitivityPatchIDs_
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

FIBase::FIBase
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleVars& primalVars,
    incompressibleAdjointVars& adjointVars,
    objectiveManager& objectiveManager,
    fv::optionAdjointList& fvOptionsAdjoint
)
:
    adjointSensitivity
    (
        mesh,
        dict,
        primalVars,
        adjointVars,
        objectiveManager,
        fvOptionsAdjoint
    ),
    shapeSensitivitiesBase(mesh, dict),
    gradDxDbMult_
    (
        IOobject
        (
            "gradDxDbMult",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(sqr(dimLength)/pow3(dimTime), Zero)
    ),
    divDxDbMult_(mesh_.nCells(), Zero),
    optionsDxDbMult_(mesh_.nCells(), Zero),
    dSfdbMult_(createZeroBoundaryPtr<vector>(mesh_)),
    dnfdbMult_(createZeroBoundaryPtr<vector>(mesh_)),
    dxdbDirectMult_(createZeroBoundaryPtr<vector>(mesh_)),

    includeDistance_(false),
    eikonalSolver_(nullptr)
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool FIBase::readDict(const dictionary& dict)
{
    if (sensitivity::readDict(dict))
    {
        if (eikonalSolver_.valid())
        {
            eikonalSolver_().readDict(dict);
        }

        return true;
    }

    return false;
}


void FIBase::accumulateIntegrand(const scalar dt)
{
    // Accumulate multiplier of grad(dxdb)
    gradDxDbMult_ += computeGradDxDbMultiplier()().T()*dt;

    // Accumulate multiplier of div(dxdb)
    PtrList<objective>& functions(objectiveManager_.getObjectiveFunctions());
    for (objective& func : functions)
    {
        divDxDbMult_ +=
            func.weight()*func.divDxDbMultiplier().primitiveField()*dt;
    }

    // Terms from fvOptions
    for (fv::optionAdjoint& optionAdj : fvOptionsAdjoint_)
    {
        optionsDxDbMult_ +=
            optionAdj.dxdbMult(adjointVars_)().primitiveField()*dt;
    }

    // Accumulate source for the adjoint to the eikonal equation
    if (includeDistance_)
    {
        eikonalSolver_->accumulateIntegrand(dt);
    }

    // Accumulate direct sensitivities
    for (const label patchI : sensitivityPatchIDs_)
    {
        const scalarField magSfDt(mesh_.boundary()[patchI].magSf()*dt);
        for (objective& func : functions)
        {
            const scalar wei = func.weight();
            dSfdbMult_()[patchI] += wei*func.dSdbMultiplier(patchI)*dt;
            dnfdbMult_()[patchI] += wei*func.dndbMultiplier(patchI)*magSfDt;
            dxdbDirectMult_()[patchI] +=
                wei*func.dxdbDirectMultiplier(patchI)*magSfDt;
        }
    }
}


void FIBase::clearSensitivities()
{
    gradDxDbMult_ = dimensionedTensor(gradDxDbMult_.dimensions(), Zero);
    divDxDbMult_ = Zero;
    optionsDxDbMult_ = vector::zero;
    dSfdbMult_() = vector::zero;
    dnfdbMult_() = vector::zero;
    dxdbDirectMult_() = vector::zero;

    if (includeDistance_)
    {
        eikonalSolver_->reset();
    }

    adjointSensitivity::clearSensitivities();
    shapeSensitivitiesBase::clear();
}


void FIBase::write(const word& baseName)
{
    adjointSensitivity::write(baseName);
    shapeSensitivitiesBase::write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

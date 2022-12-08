/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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

#include "objectiveNutSqr.H"
#include "incompressibleAdjointSolver.H"
#include "createZeroField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveNutSqr, 0);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectiveNutSqr,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveNutSqr::objectiveNutSqr
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    zones_(mesh_.cellZones().indices(dict.get<wordRes>("zones")))
{
    // Check if cellZones provided include at least one cell
    checkCellZonesSize(zones_);
    // Allocate dJdTMvar1Ptr_ and dJdTMvar2Ptr_ if needed
    allocatedJdTurbulence();
    // Allocate term to be added to volume-based sensitivity derivatives
    divDxDbMultPtr_.reset
    (
        createZeroFieldPtr<scalar>
        (
            mesh_,
            ("divDxdbMult"+type()) ,
            // Dimensions are set in a way that the gradient of this term
            // matches the source of the adjoint grid displacement PDE
            sqr(dimLength)/pow3(dimTime)
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar objectiveNutSqr::J()
{
    J_ = Zero;

    const autoPtr<incompressible::RASModelVariables>&
       turbVars = vars_.RASModelVariables();
    const volScalarField& nut = turbVars->nutRefInst();

    //scalar zoneVol(Zero);
    for (const label zI : zones_)
    {
        const cellZone& zoneI = mesh_.cellZones()[zI];
        for (const label cellI : zoneI)
        {
            J_ += sqr(nut[cellI])*(mesh_.V()[cellI]);
            //zoneVol += mesh_.V()[cellI];
        }
    }
    reduce(J_, sumOp<scalar>());
    //reduce(zoneVol, sumOp<scalar>());

    return J_;
}


void objectiveNutSqr::update_dJdv()
{
    // Add source from possible dependencies of nut on U
    if (mesh_.foundObject<incompressibleAdjointSolver>(adjointSolverName_))
    {
        const incompressibleAdjointSolver& adjSolver =
            mesh_.lookupObject<incompressibleAdjointSolver>
                (adjointSolverName_);
        const autoPtr<incompressibleAdjoint::adjointRASModel>& adjointRAS =
            adjSolver.getAdjointVars().adjointTurbulence();
        const autoPtr<incompressible::RASModelVariables>& turbVars =
            vars_.RASModelVariables();
        tmp<volScalarField> dnutdUMult = 2*turbVars->nutRef();
        tmp<volVectorField> dnutdU = adjointRAS->nutJacobianU(dnutdUMult);
        if (dnutdU)
        {
            if (!dJdvPtr_)
            {
                dJdvPtr_.reset
                (
                    createZeroFieldPtr<vector>
                    (
                        mesh_,
                        "dJdv" + type(),
                        dimLength/sqr(dimTime)
                    )
                );
            }
            for (const label zI : zones_)
            {
                const cellZone& zoneI = mesh_.cellZones()[zI];
                for (const label cellI : zoneI)
                {
                    dJdvPtr_()[cellI] = dnutdU()[cellI];
                }
            }
        }
    }
}


void objectiveNutSqr::update_dJdTMvar1()
{
    const autoPtr<incompressible::RASModelVariables>&
       turbVars = vars_.RASModelVariables();
    const volScalarField& nut = turbVars->nutRef();
    volScalarField JacobianMultiplier(2*nut);

    update_dJdTMvar
    (
        dJdTMvar1Ptr_,
        &incompressibleAdjoint::adjointRASModel::nutJacobianTMVar1,
        JacobianMultiplier,
        zones_
    );
}


void objectiveNutSqr::update_dJdTMvar2()
{
    const autoPtr<incompressible::RASModelVariables>&
       turbVars = vars_.RASModelVariables();
    const volScalarField& nut = turbVars->nutRef();
    volScalarField JacobianMultiplier(2*nut);

    update_dJdTMvar
    (
        dJdTMvar2Ptr_,
        &incompressibleAdjoint::adjointRASModel::nutJacobianTMVar2,
        JacobianMultiplier,
        zones_
    );
}


void objectiveNutSqr::update_divDxDbMultiplier()
{
    const autoPtr<incompressible::RASModelVariables>&
       turbVars = vars_.RASModelVariables();
    const volScalarField& nut = turbVars->nutRef();

    volScalarField& divDxDbMult = divDxDbMultPtr_();

    for (const label zI : zones_)
    {
        const cellZone& zoneI = mesh_.cellZones()[zI];
        for (const label cellI : zoneI)
        {
            divDxDbMult[cellI] = sqr(nut[cellI]);
        }
    }
    divDxDbMult.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //

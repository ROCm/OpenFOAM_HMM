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
#include "createZeroField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveNutSqr, 1);
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
    zones_
    (
        mesh_.cellZones().indices(this->dict().get<wordRes>("zones"))
    )
{
    // Allocate source term for the adjoint turbulence model
    dJdTMvar1Ptr_.reset
    (
        createZeroFieldPtr<scalar>
        (
            mesh_,
            "ATMSource",
            (dimless/dimTime/dimTime)
        )
    );
    // Allocate term to be added to volume-based sensitivity derivatives
    divDxDbMultPtr_.reset
    (
        createZeroFieldPtr<scalar>
        (
            mesh_,
            ("divDxdbMult"+type()) ,
            //variable dimensions!!
            //Dummy zeroes. Only the internalField will be used
            dimless
        )
    );
    // set file pointer
    //objFunctionFilePtr_ = new OFstream(objFunctionFolder_/type());
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


void objectiveNutSqr::update_dJdTMvar1()
{
    const autoPtr<incompressible::RASModelVariables>&
       turbVars = vars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = vars_.laminarTransport();
    const volScalarField& nut = turbVars->nutRef();

    tmp<volScalarField> tnutJacobian = turbVars->nutJacobianVar1(lamTransp);
    const volScalarField& nutJacobian = tnutJacobian();

    volScalarField& dJdTMvar1 = dJdTMvar1Ptr_();

    for (const label zI : zones_)
    {
        const cellZone& zoneI = mesh_.cellZones()[zI];
        for (const label cellI : zoneI)
        {
            dJdTMvar1[cellI] = 2.*nut[cellI]*nutJacobian[cellI];
        }
    }
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
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //

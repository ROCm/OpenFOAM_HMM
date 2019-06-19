/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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

#include "SpalartAllmaras.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASVariables
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmaras, 0);
addToRunTimeSelectionTable(RASModelVariables, SpalartAllmaras, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmaras::SpalartAllmaras
(
    const fvMesh& mesh,
    const solverControl& SolverControl
)
:
    RASModelVariables(mesh, SolverControl)
{
    hasTMVar1_ = true;
    TMVar1Ptr_ = mesh_.getObjectPtr<volScalarField>("nuTilda");
    TMVar1BaseName_ = "nuTilda";

    TMVar2Ptr_ =
        new volScalarField
        (
            IOobject
            (
                "dummySpalartAllmarasVar2",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );

    hasNut_ = true;
    nutPtr_ = mesh_.getObjectPtr<volScalarField>("nut");

    hasDist_ = true;
    dPtr_ = mesh_.getObjectPtr<volScalarField>("yWall");

    allocateInitValues();
    allocateMeanFields();
}


SpalartAllmaras::~SpalartAllmaras ()
{
    // nullify pointer
    TMVar1Ptr_ = nullptr;
    nutPtr_ = nullptr;
    dPtr_ = nullptr;

    // nullify pointer and delete allocated field
    delete TMVar2Ptr_;
    TMVar2Ptr_ = nullptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmaras::nutJacobianVar1
(
    const singlePhaseTransportModel& laminarTransport
) const
{
    tmp<volScalarField> tnutJacobian
    (
        new volScalarField
        (
            IOobject
            (
                "nutJacobianVar1",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    volScalarField& nutJacobian = tnutJacobian.ref();

    const volScalarField& nu = laminarTransport.nu();
    const volScalarField& nuTilda = TMVar1();
    volScalarField chi(nuTilda/nu);
    volScalarField chi3(pow3(chi));
    scalar Cv13 = pow3(7.1);
    volScalarField fv1(chi3/(chi3 + Cv13));
    volScalarField fv1ByChi2Sqr(sqr(chi/(chi3 + Cv13)));
    volScalarField Cdfv1(3.0*Cv13*fv1ByChi2Sqr);
    nutJacobian = Cdfv1*chi + fv1;

    return tnutJacobian;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASVariables
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

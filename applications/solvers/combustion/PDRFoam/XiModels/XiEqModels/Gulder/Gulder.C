/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

#include "Gulder.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{
    defineTypeNameAndDebug(Gulder, 0);
    addToRunTimeSelectionTable(XiEqModel, Gulder, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::Gulder::Gulder
(
    const dictionary& XiEqProperties,
    const hhuCombustionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(XiEqProperties, thermo, turbulence, Su),
    XiEqCoef_(readScalar(XiEqModelCoeffs_.lookup("XiEqCoef"))),
    SuMin_(0.01*Su.average()),
    uPrimeCoef_(readScalar(XiEqModelCoeffs_.lookup("uPrimeCoef")))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::XiEqModels::Gulder::~Gulder()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::Gulder::XiEq() const
{
    volScalarField up(sqrt((2.0/3.0)*turbulence_.k()));
    const volScalarField& epsilon = turbulence_.epsilon();
    const fvMesh& mesh = Su_.mesh();

    const volVectorField& U = mesh.lookupObject<volVectorField>("U");

    const volSymmTensorField& CT = mesh.lookupObject<volSymmTensorField>("CT");
    const volScalarField& Nv = mesh.lookupObject<volScalarField>("Nv");
    const volSymmTensorField& nsv = 
        mesh.lookupObject<volSymmTensorField>("nsv");

    tmp<volScalarField> tN
    (
        new volScalarField
        (
            IOobject
            (
                "tN",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", Nv.dimensions(), 0.0),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    volScalarField& N = tN();

    N.internalField() = Nv.internalField()*pow(mesh.V(), 2.0/3.0);

    tmp<volSymmTensorField> tns
    (
        new volSymmTensorField
        (
            IOobject
            (
                "tns",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor
            (
                "zero",
                nsv.dimensions(),
                pTraits<symmTensor>::zero
            )
        )
    );

    volSymmTensorField& ns = tns();

    ns.internalField() = nsv.internalField()*pow(mesh.V(), 2.0/3.0);

    const volVectorField Uhat
    (
        U/(mag(U) + dimensionedScalar("Usmall", U.dimensions(), 1e-4))
    );

    const volScalarField nr(sqrt(max(N - (Uhat & ns & Uhat), 1e-4)));

    const scalarField cellWidth(pow(mesh.V(), 1.0/3.0));

    const scalarField upLocal(uPrimeCoef_*sqrt((U & CT & U)*cellWidth));

    const scalarField deltaUp(upLocal*(max(1.0, pow(nr, 0.5)) - 1.0));

    up.internalField() += deltaUp;

    volScalarField tauEta(sqrt(mag(thermo_.muu()/(thermo_.rhou()*epsilon))));

    volScalarField Reta =
    (
        up
      /
        (
            sqrt(epsilon*tauEta)
          + dimensionedScalar("1e-8", up.dimensions(), 1e-8)
        )
    );

    return (1.0 + XiEqCoef_*sqrt(up/(Su_ + SuMin_))*Reta);
}


bool Foam::XiEqModels::Gulder::read(const dictionary& XiEqProperties)
{
    XiEqModel::read(XiEqProperties);

    XiEqModelCoeffs_.lookup("XiEqCoef") >> XiEqCoef_;
    XiEqModelCoeffs_.lookup("uPrimeCoef") >> uPrimeCoef_;

    return true;
}


// ************************************************************************* //

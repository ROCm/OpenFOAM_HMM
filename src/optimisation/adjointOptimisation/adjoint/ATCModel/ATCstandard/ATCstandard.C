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

#include "ATCstandard.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ATCstandard, 0);
addToRunTimeSelectionTable
(
    ATCModel,
    ATCstandard,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ATCstandard::ATCstandard
(
    const fvMesh& mesh,
    const incompressibleVars& primalVars,
    const incompressibleAdjointVars& adjointVars,
    const dictionary& dict
)
:
    ATCModel(mesh, primalVars, adjointVars, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ATCstandard::addATC(fvVectorMatrix& UaEqn)
{
    const volVectorField& U = primalVars_.U();
    const volVectorField& Ua = adjointVars_.UaInst();
    const surfaceScalarField& phi = primalVars_.phi();

    // Build U to go into the ATC term, based on whether to smooth field or not
    autoPtr<volVectorField> UForATC(nullptr);
    if (reconstructGradients_)
    {
        UForATC.reset(new volVectorField(fvc::reconstruct(phi)));
    }
    else
    {
        UForATC.reset(new volVectorField(U));
    }

    // Main ATC term
    ATC_ = (fvc::grad(UForATC(), "gradUATC") & Ua);

    if (extraConvection_ > 0)
    {
        // Implicit part added to increase diagonal dominance
        // Note: Maybe this needs to be multiplied with the ATClimiter ??
        UaEqn += extraConvection_*fvm::div(-phi, Ua);

        // correct rhs due to implicitly augmenting the adjoint convection
        ATC_ += extraConvection_*(fvc::grad(Ua, "gradUaATC")().T() & U);
    }

    //zero ATC on cells next to given patch types
    smoothATC();

    // actual ATC term
    UaEqn += fvm::Su(ATC_, Ua);
}


tmp<volTensorField> ATCstandard::getFISensitivityTerm() const
{
    tmp<volTensorField> tvolSDTerm
    (
        new volTensorField
        (
            IOobject
            (
                "ATCFISensitivityTerm" + type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor(sqr(dimLength)/pow(dimTime, 3), Zero)
        )
    );

    volTensorField& volSDTerm = tvolSDTerm.ref();

    const volVectorField& U = primalVars_.U();
    const volVectorField& Ua = adjointVars_.Ua();

    volTensorField gradU(fvc::grad(U));

    // Explicitly correct the boundary gradient to get rid of the
    // tangential component
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();
            const vectorField& nf = tnf();
            gradU.boundaryFieldRef()[patchI] =
                nf*U.boundaryField()[patchI].snGrad();
        }
    }

    const volScalarField& mask = getLimiter();

    volSDTerm = -(gradU & Ua)*U*mask;

    return tvolSDTerm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "LRRDiffStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESmodels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LRRDiffStress, 0);
addToRunTimeSelectionTable(LESmodel, LRRDiffStress, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
LRRDiffStress::LRRDiffStress
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESmodel(typeName, U, phi, transport),
    GenSGSStress(U, phi, transport),

    ck_(LESmodelProperties().lookupOrAddDefault<scalar>("ck", 0.09)),
    c1_(LESmodelProperties().lookupOrAddDefault<scalar>("c1", 1.8)),
    c2_(LESmodelProperties().lookupOrAddDefault<scalar>("c2", 0.6))
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void LRRDiffStress::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();

    GenSGSStress::correct(gradU);

    volSymmTensorField D = symm(gradU);

    volSymmTensorField P = -twoSymm(B_ & gradU);

    volScalarField K = 0.5*tr(B_);
    volScalarField Epsilon = 2*nuEff()*magSqr(D);

    solve
    (
        fvm::ddt(B_)
      + fvm::div(phi(), B_)
      - fvm::laplacian(DBEff(), B_)
      + fvm::Sp(c1_*Epsilon/K, B_)
     ==
        P
      - (0.667*(1.0 - c1_)*I)*Epsilon
      - c2_*(P - 0.333*I*tr(P))
      - (0.667 - 2*c1_)*I*pow(K, 1.5)/delta()
    );


    // Bounding the component kinetic energies

    forAll(B_, celli)
    {
        B_[celli].component(tensor::XX) =
            max(B_[celli].component(tensor::XX), k0().value());
        B_[celli].component(tensor::YY) =
            max(B_[celli].component(tensor::YY), k0().value());
        B_[celli].component(tensor::ZZ) =
            max(B_[celli].component(tensor::ZZ), k0().value());
    }

    K = 0.5*tr(B_);
    bound(K, k0());

    nuSgs_ = ck_*sqrt(K)*delta();
    nuSgs_.correctBoundaryConditions();
}


bool LRRDiffStress::read()
{
    if (GenSGSStress::read())
    {
        LESmodelProperties().readIfPresent<scalar>("ck", ck_);
        LESmodelProperties().readIfPresent<scalar>("c1", c1_);
        LESmodelProperties().readIfPresent<scalar>("c2", c2_);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESmodels
} // End namespace Foam

// ************************************************************************* //

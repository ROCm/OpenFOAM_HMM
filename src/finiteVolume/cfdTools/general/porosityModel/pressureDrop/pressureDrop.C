/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "pressureDrop.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(pressureDrop, 0);
        addToRunTimeSelectionTable(porosityModel, pressureDrop, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::pressureDrop::pressureDrop
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),
    mDotvsDp_(DataEntry<scalar>::New("mDotvsDp", coeffs_)),
    lRef_(readScalar(coeffs_.lookup("lRef"))),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::pressureDrop::~pressureDrop()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::pressureDrop::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    scalar rhoScale = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        apply(Udiag, Usource, V, rho, U, rhoScale);
    }
    else
    {
        coeffs_.lookup("rhoRef") >> rhoScale;
        apply(Udiag, Usource, V, geometricOneField(), U, rhoScale);
    }

}


void Foam::porosityModels::pressureDrop::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField&
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    apply(Udiag, Usource, V, rho, U, 1.0);
}


void Foam::porosityModels::pressureDrop::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();

    scalar rhoScale = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        apply(AU, rho, U, rhoScale);
    }
    else
    {
        coeffs_.lookup("rhoRef") >> rhoScale;
        apply(AU, geometricOneField(), U, rhoScale);
    }
}


void Foam::porosityModels::pressureDrop::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


// ************************************************************************* //

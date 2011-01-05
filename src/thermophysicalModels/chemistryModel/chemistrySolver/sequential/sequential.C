/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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

#include "sequential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::sequential<CompType, ThermoType>::sequential
(
    ODEChemistryModel<CompType, ThermoType>& model,
    const word& modelName
)
:
    chemistrySolver<CompType, ThermoType>(model, modelName),
    coeffsDict_(model.subDict(modelName + "Coeffs")),
    cTauChem_(readScalar(coeffsDict_.lookup("cTauChem"))),
    eqRateLimiter_(coeffsDict_.lookup("equilibriumRateLimiter"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::sequential<CompType, ThermoType>::~sequential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::scalar Foam::sequential<CompType, ThermoType>::solve
(
    scalarField &c,
    const scalar T,
    const scalar p,
    const scalar t0,
    const scalar dt
) const
{
    scalar tChemInv = SMALL;

    scalar pf, cf, pb, cb;
    label lRef, rRef;

    forAll(this->model_.reactions(), i)
    {
        const Reaction<ThermoType>& R = this->model_.reactions()[i];

        scalar omega = this->model_.omega
        (
            R, c, T, p, pf, cf, lRef, pb, cb, rRef
        );

        if (eqRateLimiter_)
        {
            if (omega < 0.0)
            {
                omega /= 1.0 + pb*dt;
            }
            else
            {
                omega /= 1.0 + pf*dt;
            }
        }

        tChemInv = max(tChemInv, mag(omega));


        // update species
        forAll(R.lhs(), specieI)
        {
            const label id = R.lhs()[specieI].index;
            const scalar sc = R.lhs()[specieI].stoichCoeff;
            c[id] -= dt*sc*omega;
            c[id] = max(0.0, c[id]);
        }

        forAll(R.rhs(), specieI)
        {
            const label id = R.rhs()[specieI].index;
            const scalar sc = R.rhs()[specieI].stoichCoeff;
            c[id] += dt*sc*omega;
            c[id] = max(0.0, c[id]);
        }
    }

    return cTauChem_/tChemInv;
}


// ************************************************************************* //

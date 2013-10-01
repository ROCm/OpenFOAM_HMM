/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "EulerImplicit.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::EulerImplicit<ChemistryModel>::EulerImplicit
(
    const fvMesh& mesh
)
:
    chemistrySolver<ChemistryModel>(mesh),
    coeffsDict_(this->subDict("EulerImplicitCoeffs")),
    cTauChem_(readScalar(coeffsDict_.lookup("cTauChem"))),
    eqRateLimiter_(coeffsDict_.lookup("equilibriumRateLimiter")),
    cTp_(this->nEqns())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::EulerImplicit<ChemistryModel>::~EulerImplicit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::EulerImplicit<ChemistryModel>::solve
(
    scalarField& c,
    scalar& T,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    deltaT = min(deltaT, subDeltaT);

    const label nSpecie = this->nSpecie();
    simpleMatrix<scalar> RR(nSpecie, 0, 0);

    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
    }

    // Calculate the absolute enthalpy
    scalar cTot = sum(c);
    typename ChemistryModel::thermoType mixture
    (
        (c[0]/cTot)*this->specieThermo_[0]
    );
    for (label i=1; i<nSpecie; i++)
    {
        mixture += (c[i]/cTot)*this->specieThermo_[i];
    }
    scalar ha = mixture.Ha(p, T);

    for (label i=0; i<nSpecie; i++)
    {
        RR.source()[i] = c[i]/deltaT;
    }

    forAll(this->reactions(), i)
    {
        scalar pf, cf, pr, cr;
        label lRef, rRef;

        scalar omegai = this->omegaI(i, c, T, p, pf, cf, lRef, pr, cr, rRef);

        scalar corr = 1.0;
        if (eqRateLimiter_)
        {
            if (omegai < 0.0)
            {
                corr = 1.0/(1.0 + pr*deltaT);
            }
            else
            {
                corr = 1.0/(1.0 + pf*deltaT);
            }
        }

        this->updateRRInReactionI(i, pr, pf, corr, lRef, rRef, p, T, RR);
    }


    for (label i=0; i<nSpecie; i++)
    {
        RR[i][i] += 1.0/deltaT;
    }

    c = RR.LUsolve();
    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
    }

    // Update the temperature
    cTot = sum(c);
    mixture = (c[0]/cTot)*this->specieThermo_[0];
    for (label i=1; i<nSpecie; i++)
    {
        mixture += (c[i]/cTot)*this->specieThermo_[i];
    }
    T = mixture.THa(ha, p, T);

    // Estimate the next time step
    scalar tMin = GREAT;
    const label nEqns = this->nEqns();

    for (label i=0; i<nSpecie; i++)
    {
        cTp_[i] = c[i];
    }
    cTp_[nSpecie] = T;
    cTp_[nSpecie+1] = p;

    scalarField dcdt(nEqns, 0.0);
    this->derivatives(0.0, cTp_, dcdt);

    const scalar sumC = sum(c);

    for (label i=0; i<nSpecie; i++)
    {
        scalar d = dcdt[i];

        if (d < -SMALL)
        {
            tMin = min(tMin, -(c[i] + SMALL)/d);
        }
        else
        {
            d = max(d, SMALL);
            scalar cm = max(sumC - c[i], 1.0e-5);
            tMin = min(tMin, cm/d);
        }
    }

    subDeltaT = cTauChem_*tMin;
}


// ************************************************************************* //

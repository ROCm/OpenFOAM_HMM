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

#include "sequential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::sequential<ChemistryModel>::sequential
(
    const fvMesh& mesh
)
:
    chemistrySolver<ChemistryModel>(mesh),
    coeffsDict_(this->subDict("sequentialCoeffs")),
    cTauChem_(readScalar(coeffsDict_.lookup("cTauChem"))),
    eqRateLimiter_(coeffsDict_.lookup("equilibriumRateLimiter"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::sequential<ChemistryModel>::~sequential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::sequential<ChemistryModel>::solve
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

    scalar tChemInv = SMALL;

    forAll(this->reactions(), i)
    {
        scalar pf, cf, pb, cb;
        label lRef, rRef;

        scalar omega = this->omegaI(i, c, T, p, pf, cf, lRef, pb, cb, rRef);

        if (eqRateLimiter_)
        {
            if (omega < 0.0)
            {
                omega /= 1.0 + pb*deltaT;
            }
            else
            {
                omega /= 1.0 + pf*deltaT;
            }
        }

        tChemInv = max(tChemInv, mag(omega));

        this->updateConcsInReactionI(i, deltaT, omega, p, T, c);
    }

    // Update the temperature
    cTot = sum(c);
    mixture = (c[0]/cTot)*this->specieThermo_[0];
    for (label i=1; i<nSpecie; i++)
    {
        mixture += (c[i]/cTot)*this->specieThermo_[i];
    }
    T = mixture.THa(ha, p, T);

    subDeltaT = cTauChem_/tChemInv;
}


// ************************************************************************* //

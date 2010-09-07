/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "standardPhaseChange.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoSingleLayer.H"
#include "specie.H"
#include "heatTransferModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(standardPhaseChange, 0);
        addToRunTimeSelectionTable
        (
            phaseChangeModel,
            standardPhaseChange,
            dictionary
        );
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::surfaceFilmModels::standardPhaseChange::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    if (Re < 5.0E+05)
    {
        return 0.664*sqrt(Re)*cbrt(Sc);
    }
    else
    {
        return 0.037*pow(Re, 0.8)*cbrt(Sc);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::standardPhaseChange::standardPhaseChange
(
    const surfaceFilmModel& owner,
    const dictionary& dict
)
:
    phaseChangeModel(typeName, owner, dict),
    deltaMin_(readScalar(coeffs_.lookup("deltaMin"))),
    L_(readScalar(coeffs_.lookup("L")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::standardPhaseChange::~standardPhaseChange()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceFilmModels::standardPhaseChange::correct
(
    const scalar dt,
    scalarField& dMass
)
{
    const thermoSingleLayer& film = refCast<const thermoSingleLayer>(owner_);

    // set local thermo properties
    const SLGThermo& thermo = film.thermo();
    const label liqId = film.liquidId();
    const liquid& liq = thermo.liquids().properties()[liqId];
    const label vapId = thermo.carrierId(thermo.liquids().components()[liqId]);

    // retrieve fields from film model
    const scalarField& delta = film.delta();
    const scalarField& YInf = film.YPrimary()[vapId];
    const scalarField& pInf = film.pPrimary();
    const scalarField& T = film.T();
    const scalarField& Ts = film.Ts();
    const scalarField& Tw = film.Tw();
    const scalarField& TInf = film.TPrimary();
    const scalarField& rho = film.rho();
    const scalarField& mu = film.mu();
    const scalarField& hs = film.hs();
    const scalarField& magSf = film.magSf();
    const scalarField hInf = film.htcs().h();
    const scalarField hFilm = film.htcw().h();
    const vectorField dU = film.UPrimary() - film.Us();

    // Reynolds number
    const scalarField Re = rho*mag(dU)*L_/mu;

    // molecular weight of vapour
    const scalar Wvap = thermo.carrier().W(vapId);

    // molecular weight of liquid
    const scalar Wliq = liq.W();

    forAll(dMass, cellI)
    {
        if (delta[cellI] > deltaMin_)
        {
            // cell pressure
            const scalar pc = pInf[cellI];

            // saturation pressure
            const scalar pSat = liq.pv(pc, Ts[cellI]);

            // calculate mass transfer
            if (pSat > pc)
            {
                // boiling
                const scalar qDotInf = hInf[cellI]*(TInf[cellI] - T[cellI]);
                const scalar qDotFilm = hFilm[cellI]*(T[cellI] - Tw[cellI]);
                dMass[cellI] +=
                    max
                    (
                        0.0,
                        dt*magSf[cellI]/hs[cellI]*(qDotInf - qDotFilm)
                    );
            }
            else
            {
                // vapour mass fraction at interface
                const scalar Ys = Wliq*pSat/(Wliq*pSat + Wvap*(pc - pSat));

                // bulk gas average density
                const scalar rhoAve = pc/(specie::RR*Ts[cellI]);

                // vapour diffusivity [m2/s]
                const scalar Dab = liq.D(pc, Ts[cellI]);

                // Schmidt number
                const scalar Sc = mu[cellI]/(rho[cellI]*(Dab + ROOTVSMALL));

                // Sherwood number
                const scalar Sh = this->Sh(Re[cellI], Sc);

                // mass transfer coefficient [m/s]
                const scalar hm = Sh*Dab/L_;

                // add mass contribution to source
                dMass[cellI] =
                    max
                    (
                        0.0,
                        dt*magSf[cellI]*rhoAve*hm*(Ys - YInf[cellI])/(1.0 - Ys)
                    );
            }
        }
    }
}


// ************************************************************************* //

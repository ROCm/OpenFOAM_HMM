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
    Tb_(readScalar(coeffs_.lookup("Tb"))),
    deltaMin_(readScalar(coeffs_.lookup("deltaMin"))),
    L_(readScalar(coeffs_.lookup("L"))),
    totalMass_(0.0),
    vapourRate_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::standardPhaseChange::~standardPhaseChange()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceFilmModels::standardPhaseChange::correct
(
    const scalar dt,
    scalarField& dMass,
    scalarField& dEnergy
)
{
    dMass = 0.0;
    dEnergy = 0.0;

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
    const scalarField& Tw = film.Tw();
    const scalarField& rho = film.rho();
    const scalarField& TInf = film.TPrimary();
    const scalarField& rhoInf = film.rhoPrimary();
    const scalarField& muInf = film.muPrimary();
    const scalarField& magSf = film.magSf();
    const scalarField hInf = film.htcs().h();
    const scalarField hFilm = film.htcw().h();
    const vectorField dU = film.UPrimary() - film.Us();
    const scalarField availableMass = (delta - deltaMin_)*rho*magSf;


    forAll(dMass, cellI)
    {
        if (delta[cellI] > deltaMin_)
        {
            // cell pressure [Pa]
            const scalar pc = pInf[cellI];

            // local temperature - impose lower limit of 200 K for stability
            const scalar Tloc = min(liq.Tc(), max(200.0, T[cellI]));

            // saturation pressure [Pa]
            const scalar pSat = liq.pv(pc, Tloc);

            // latent heat [J/kg]
            const scalar hVap = liq.hl(pc, Tloc);

            // calculate mass transfer
            if (pSat >= 0.95*pc)
            {
                // boiling
                const scalar qDotInf = hInf[cellI]*(TInf[cellI] - T[cellI]);
                const scalar qDotFilm = hFilm[cellI]*(T[cellI] - Tw[cellI]);

                const scalar Cp = liq.Cp(pc, Tloc);
                const scalar qCorr = availableMass[cellI]*Cp*(T[cellI] - Tb_);
                dMass[cellI] =
                    dt*magSf[cellI]/hVap*(qDotInf + qDotFilm)
                  + qCorr/hVap;

            }
            else
            {
                // Primary region density [kg/m3]
                const scalar rhoInfc = rhoInf[cellI];

                // Primary region viscosity [Pa.s]
                const scalar muInfc = muInf[cellI];

                // Reynolds number
                const scalar Re = rhoInfc*mag(dU[cellI])*L_/muInfc;

                // molecular weight of vapour [kg/kmol]
                const scalar Wvap = thermo.carrier().W(vapId);

                // molecular weight of liquid [kg/kmol]
                const scalar Wliq = liq.W();

                // vapour mass fraction at interface
                const scalar Ys = Wliq*pSat/(Wliq*pSat + Wvap*(pc - pSat));

                // vapour diffusivity [m2/s]
                const scalar Dab = liq.D(pc, Tloc);

                // Schmidt number
                const scalar Sc = muInfc/(rhoInfc*(Dab + ROOTVSMALL));

                // Sherwood number
                const scalar Sh = this->Sh(Re, Sc);

                // mass transfer coefficient [m/s]
                const scalar hm = Sh*Dab/L_;

                // add mass contribution to source
                dMass[cellI] =
                    dt*magSf[cellI]*rhoInfc*hm*(Ys - YInf[cellI])/(1.0 - Ys);
            }

            dMass[cellI] = min(availableMass[cellI], max(0.0, dMass[cellI]));
            dEnergy[cellI] = dMass[cellI]*hVap;
        }
    }

    const scalar sumdMass = sum(dMass);
    totalMass_ += sumdMass;
    vapourRate_ = sumdMass/owner().time().deltaTValue();
}


void Foam::surfaceFilmModels::standardPhaseChange::info() const
{
    Info<< indent << "mass phase change  = "
        << returnReduce(totalMass_, sumOp<scalar>()) << nl
        << indent << "vapourisation rate = "
        << returnReduce(vapourRate_, sumOp<scalar>()) << nl;
}


// ************************************************************************* //

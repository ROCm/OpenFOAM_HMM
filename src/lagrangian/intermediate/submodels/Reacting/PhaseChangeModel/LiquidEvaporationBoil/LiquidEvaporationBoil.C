/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "LiquidEvaporationBoil.H"
#include "specie.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::tmp<Foam::scalarField> Foam::LiquidEvaporationBoil<CloudType>::calcXc
(
    const label cellI
) const
{
    scalarField Xc(this->owner().thermo().carrier().Y().size());

    forAll(Xc, i)
    {
        Xc[i] =
            this->owner().thermo().carrier().Y()[i][cellI]
           /this->owner().thermo().carrier().W(i);
    }

    return Xc/sum(Xc);
}


template<class CloudType>
Foam::scalar Foam::LiquidEvaporationBoil<CloudType>::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    return 2.0 + 0.6*Foam::sqrt(Re)*cbrt(Sc);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LiquidEvaporationBoil<CloudType>::LiquidEvaporationBoil
(
    const dictionary& dict,
    CloudType& owner
)
:
    PhaseChangeModel<CloudType>(dict, owner, typeName),
    liquids_(owner.thermo().liquids()),
    activeLiquids_(this->coeffDict().lookup("activeLiquids")),
    liqToCarrierMap_(activeLiquids_.size(), -1),
    liqToLiqMap_(activeLiquids_.size(), -1)
{
    if (activeLiquids_.size() == 0)
    {
        WarningIn
        (
            "Foam::LiquidEvaporationBoil<CloudType>::LiquidEvaporationBoil"
            "("
                "const dictionary& dict, "
                "CloudType& owner"
            ")"
        )   << "Evaporation model selected, but no active liquids defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating liquid species:" << endl;

        // Determine mapping between liquid and carrier phase species
        forAll(activeLiquids_, i)
        {
            Info<< "    " << activeLiquids_[i] << endl;
            liqToCarrierMap_[i] =
                owner.composition().globalCarrierId(activeLiquids_[i]);
        }

        // Determine mapping between model active liquids and global liquids
        const label idLiquid = owner.composition().idLiquid();
        forAll(activeLiquids_, i)
        {
            liqToLiqMap_[i] =
                owner.composition().localId(idLiquid, activeLiquids_[i]);
        }
    }
}


template<class CloudType>
Foam::LiquidEvaporationBoil<CloudType>::LiquidEvaporationBoil
(
    const LiquidEvaporationBoil<CloudType>& pcm
)
:
    PhaseChangeModel<CloudType>(pcm),
    liquids_(pcm.owner().thermo().liquids()),
    activeLiquids_(pcm.activeLiquids_),
    liqToCarrierMap_(pcm.liqToCarrierMap_),
    liqToLiqMap_(pcm.liqToLiqMap_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LiquidEvaporationBoil<CloudType>::~LiquidEvaporationBoil()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::LiquidEvaporationBoil<CloudType>::calculate
(
    const scalar dt,
    const label cellI,
    const scalar Re,
    const scalar Pr,
    const scalar d,
    const scalar nu,
    const scalar T,
    const scalar Ts,
    const scalar pc,
    const scalar Tc,
    const scalarField& Yl,
    scalarField& dMassPC
) const
{
    // construct carrier phase species volume fractions for cell, cellI
    const scalarField XcMix(calcXc(cellI));

    // liquid volume fraction
    const scalarField X(liquids_.X(Yl));

    // droplet surface pressure
    scalar ps = liquids_.pv(pc, T, X);

    // vapour density at droplet surface [kg/m3]
    scalar rhos = ps*liquids_.W(X)/(specie::RR*Ts);

    // thermal conductivity of carrier [W/m/K]
    scalar kappac = 0.0;
    forAll(this->owner().thermo().carrier().Y(), i)
    {
        const scalar Yc = this->owner().thermo().carrier().Y()[i][cellI];
        kappac += Yc*this->owner().thermo().carrier().kappa(i, Ts);
    }

    // calculate mass transfer of each specie in liquid
    forAll(activeLiquids_, i)
    {
        const label gid = liqToCarrierMap_[i];
        const label lid = liqToLiqMap_[i];

        // boiling temperature at cell pressure for liquid species lid [K]
        const scalar TBoil = liquids_.properties()[lid].pvInvert(pc);

        // limit droplet temperature to boiling/critical temperature
        const scalar Td = min(T, 0.999*TBoil);

        // saturation pressure for liquid species lid [Pa]
        const scalar pSat = liquids_.properties()[lid].pv(pc, Td);

        // carrier phase concentration
        const scalar Xc = XcMix[gid];

        if (Xc*pc > pSat)
        {
            // saturated vapour - no phase change
        }
        else
        {
            if (pSat > 0.999*pc)
            {
                // boiling

//                const scalar deltaT = max(Tc - Td, 0.5);
                const scalar deltaT = max(Tc - T, 0.5);

                // liquid specific heat capacity
                const scalar Cp = liquids_.properties()[lid].Cp(pc, Td);

                // vapour heat of fomation
                const scalar hv = liquids_.properties()[lid].hl(pc, Td);

                // Nusselt number
                const scalar Nu = 2.0 + 0.6*sqrt(Re)*cbrt(Pr);

                const scalar lg = log(1.0 + Cp*deltaT/max(SMALL, hv));

                // mass transfer [kg]
                dMassPC[lid] += pi*kappac*Nu*lg*d/Cp*dt;
            }
            else
            {
                // evaporation

                // surface molar fraction - Raoult's Law
                const scalar Xs = X[lid]*pSat/pc;

                // molar ratio
                const scalar Xr = (Xs - Xc)/max(SMALL, 1.0 - Xs);


                if (Xr > 0)
                {
                    // vapour diffusivity [m2/s]
                    const scalar Dab = liquids_.properties()[lid].D(pc, Td);

                    // Schmidt number
                    const scalar Sc = nu/(Dab + ROOTVSMALL);

                    // Sherwood number
                    const scalar Sh = this->Sh(Re, Sc);

                    // mass transfer coefficient [m/s]
                    const scalar kc = Sh*Dab/(d + ROOTVSMALL);

                    // mass transfer [kg]
                    dMassPC[lid] += pi*sqr(d)*kc*rhos*log(1.0 + Xr)*dt;
                }
            }
        }
    }
}


template<class CloudType>
Foam::scalar Foam::LiquidEvaporationBoil<CloudType>::dh
(
    const label idc,
    const label idl,
    const label p,
    const label T
) const
{
    scalar dh = 0;

    typedef PhaseChangeModel<CloudType> parent;
    switch (parent::enthalpyTransfer_)
    {
        case (parent::etLatentHeat):
        {
            dh = liquids_.properties()[idl].hl(p, T);
            break;
        }
        case (parent::etEnthalpyDifference):
        {
            scalar hc = this->owner().composition().carrier().H(idc, T);
            scalar hp = liquids_.properties()[idl].h(p, T);

            dh = hc - hp;
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::LiquidEvaporationBoil<CloudType>::dh"
                "("
                    "const label, "
                    "const label, "
                    "const label, "
                    "const label"
                ")"
            )   << "Unknown enthalpyTransfer type" << abort(FatalError);
        }
    }

    return dh;
}


// ************************************************************************* //

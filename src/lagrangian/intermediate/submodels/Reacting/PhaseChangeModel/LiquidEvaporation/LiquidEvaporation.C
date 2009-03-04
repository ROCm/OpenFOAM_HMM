/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "LiquidEvaporation.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class CloudType>
Foam::label Foam::LiquidEvaporation<CloudType>::carrierSpecieId
(
    const word& specieName
) const
{
    forAll (this->owner().carrierThermo().composition().Y(), i)
    {
        if
        (
            this->owner().carrierThermo().composition().Y()[i].name()
         == specieName
        )
        {
            return i;
        }
    }

    wordList species(this->owner().carrierThermo().composition().Y().size());
    forAll (this->owner().carrierThermo().composition().Y(), i)
    {
        species[i] = this->owner().carrierThermo().composition().Y()[i].name();
    }

    FatalErrorIn
    (
        "Foam::label Foam::LiquidEvaporation<CloudType>::carrierSpecieId"
        "("
            "const word&"
        ") const"
    )   << "Could not find " << specieName << " in species list" << nl
        <<  "Avialable species:" << nl << species
        << nl << exit(FatalError);

    return -1;
}


template <class CloudType>
Foam::scalar Foam::LiquidEvaporation<CloudType>::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    return 2.0 + 0.6*Foam::sqrt(Re)*pow(Sc, 0.333);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LiquidEvaporation<CloudType>::LiquidEvaporation
(
    const dictionary& dict,
    CloudType& owner
)
:
    PhaseChangeModel<CloudType>(dict, owner, typeName),
    liquids_
    (
        liquidMixture::New
        (
            owner.mesh().objectRegistry::lookupObject<dictionary>
            (
                owner.carrierThermo().name()
            )
        )
    ),
    Tvap_(readScalar(this->coeffDict().lookup("Tvap"))),
    evapProps_(this->coeffDict().lookup("activeLiquids")),
    liqToGasMap_(evapProps_.size())
{
    if (evapProps_.size() == 0)
    {
        WarningIn
        (
            "Foam::LiquidEvaporation<CloudType>::LiquidEvaporation"
            "("
                "const dictionary& dict, "
                "CloudType& owner"
            ")"
        )   << "Evaporation model selected, but no active liquids defined"
            << nl << endl;
    }

    // Calculate mapping between liquid and carrier phase species
    forAll(evapProps_, i)
    {
        liqToGasMap_[i] = carrierSpecieId(evapProps_[i].name());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LiquidEvaporation<CloudType>::~LiquidEvaporation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::LiquidEvaporation<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::LiquidEvaporation<CloudType>::calculate
(
    const scalar T,
    const scalar d,
    const scalarField& Xc,
    scalarList& dMassMT,
    const vector& Ur,
    const scalar Tc,
    const scalar pc,
    const scalar nuc,
    const scalar dt
) const
{
    scalar dMassTot = 0.0;

    if (T < Tvap_)
    {
        // not reached model activation temperature
        return dMassTot;
    }
    else
    {
        // droplet area
        scalar A = mathematicalConstant::pi*sqr(d);

        // Reynolds number
        scalar Re = mag(Ur)*d/(nuc + ROOTVSMALL);

        // Calculate mass transfer of each specie in liquid
        forAll(evapProps_, i)
        {
            // Diffusion coefficient for species i
            scalar Dab = evapProps_[i].Dab();

            // Saturation pressure for species i at temperature T
            scalar pSat = evapProps_[i].TvsPSat().value(T);

            // Schmidt number
            scalar Sc = nuc/(Dab + ROOTVSMALL);

            // Sherwood number
            scalar Sh = this->Sh(Re, Sc);

            // mass transfer coefficient [m/s]
            scalar kc = Sh*Dab/(d + ROOTVSMALL);

            // vapour concentration at droplet surface [kgmol/m3]
            scalar Cs = pSat/(specie::RR*T);

            // vapour concentration in bulk gas [kgmol/m3]
            scalar Cinf = Xc[i]*pc/(specie::RR*Tc);

            // molar flux of vapour [kgmol/m2/s]
            scalar Ni = max(kc*(Cs - Cinf), 0.0);

            // mass transfer
            label globalLiqId = liqToGasMap_[i];
            scalar dm = Ni*A*liquids_->properties()[globalLiqId].W()*dt;
            dMassMT[globalLiqId] -= dm;
            dMassTot += dm;
        }
    }

    return dMassTot;
}


// ************************************************************************* //

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


template<class CloudType>
Foam::scalarField Foam::LiquidEvaporation<CloudType>::calcXc
(
    const label cellI
) const
{
    scalarField Xc(this->owner().carrierThermo().composition().Y().size());

    scalar Winv = 0.0;
    forAll(Xc, i)
    {
        scalar Y = this->owner().carrierThermo().composition().Y()[i][cellI];
        Winv += Y/this->owner().gases()[i].W();
        Xc[i] = Y/this->owner().gases()[i].W();
    }

    return Xc/Winv;
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
    activeLiquids_(this->coeffDict().lookup("activeLiquids")),
    liqToGasMap_(activeLiquids_.size(), -1),
    liqToLiqMap_(activeLiquids_.size(), -1)
{
    if (activeLiquids_.size() == 0)
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

    // Determine mapping between liquid and carrier phase species
    forAll(activeLiquids_, i)
    {
        liqToGasMap_[i] = carrierSpecieId(activeLiquids_[i]);
    }

    // Determine mapping between local and global liquids
    forAll(activeLiquids_, i)
    {
        forAll(liquids_->components(), j)
        {
            if (liquids_->components()[j] == activeLiquids_[i])
            {
                liqToLiqMap_[i] = j;
            }
        }

        if (liqToLiqMap_[i] < 0)
        {
            FatalErrorIn
            (
                "Foam::LiquidEvaporation<CloudType>::LiquidEvaporation"
                "("
                    "const dictionary& dict, "
                    "CloudType& owner"
                ")"
            )   << "Unable to find liquid species " << activeLiquids_[i]
                << " in global liquids list. Avaliable liquids:" << nl
                << liquids_->components() << nl << exit(FatalError);
        }
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
    const scalar dt,
    const label cellI,
    const scalar T,
    const scalar pc,
    const scalar d,
    const scalar Tc,
    const scalar nuc,
    const vector& Ur,
    scalarList& dMassMT
) const
{
    // initialise total mass transferred from the particle to carrier phase
    scalar dMassTot = 0.0;

    // construct carrier phase species volume fractions for cell, cellI
    scalarField Xc = calcXc(cellI);

    // droplet surface area
    scalar A = mathematicalConstant::pi*sqr(d);

    // Reynolds number
    scalar Re = mag(Ur)*d/(nuc + ROOTVSMALL);

    // calculate mass transfer of each specie in liquid
    forAll(activeLiquids_, i)
    {
        label gid = liqToGasMap_[i];
        label lid = liqToLiqMap_[i];

        // vapour diffusivity [m2/s]
        scalar Dab = liquids_->properties()[lid].D(pc, T);

        // saturation pressure for species i [pa]
        scalar pSat = liquids_->properties()[lid].pv(pc, T);

        // Schmidt number
        scalar Sc = nuc/(Dab + ROOTVSMALL);

        // Sherwood number
        scalar Sh = this->Sh(Re, Sc);

        // mass transfer coefficient [m/s]
        scalar kc = Sh*Dab/(d + ROOTVSMALL);

        // vapour concentration at droplet surface [kgmol/m3]
        scalar Cs = pSat/(specie::RR*T);

        // vapour concentration in bulk gas [kgmol/m3]
        scalar Cinf = Xc[gid]*pc/(specie::RR*Tc);

        // molar flux of vapour [kgmol/m2/s]
        scalar Ni = max(kc*(Cs - Cinf), 0.0);

        // mass transfer [kg]
        scalar dm = Ni*A*liquids_->properties()[lid].W()*dt;
        dMassMT[gid] += dm;
        dMassTot += dm;
    }

    return dMassTot;
}


// ************************************************************************* //

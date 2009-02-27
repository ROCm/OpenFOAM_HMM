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

#include "liquidEvaporation.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class CloudType>
Foam::scalar Foam::liquidEvaporation<CloudType>::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    return 2.0 + 0.6*Foam::sqrt(Re)*pow(Sc, 0.333);
}


template <class CloudType>
Foam::scalar Foam::liquidEvaporation<CloudType>::pSat
(
    const label i,
    const scalar T
) const
{
    const List<Tuple2<scalar, scalar> >& pSat = pSat_[i];

    label id = 0;
    label nT = pSat.size();
    while ((id < nT) && (pSat[id].first() < T))
    {
        id++;
    }

    if (id == 0 || id == nT - 1)
    {
        return pSat[id].second();
    }
    else
    {
        return
            (pSat[id].first() - T)
           /(pSat[id].first() - pSat[id-1].first())
           *(pSat[id].second() - pSat[id-1].second())
          + pSat[id-1].second();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::liquidEvaporation<CloudType>::liquidEvaporation
(
    const dictionary& dict,
    CloudType& owner
)
:
    PhaseChangeModel<CloudType>(dict, owner, typeName),
    gases_(owner.gases()),
    liquids_
    (
        liquidMixture::New
        (
            owner.mesh().objectRegistry::lookupObject<dictionary>
            (
                owner.carrierThermo().name().name()
            )
        )
    ),
    Tvap_(readScalar(coeffDict().lookup("Tvap"))),
    Dm_(coeffDict().lookup("DiffusionCoeffs")),
    pSat_(coeffDict().lookup("pSat"))
{
    if (liquids_.size() != Dm_.size())
    {
        FatalErrorIn
        (
            "Foam::liquidEvaporation<CloudType>::liquidEvaporation\n"
            "(\n"
            "    const dictionary& dict,\n"
            "    CloudType& cloud\n"
            ")\n"
        )   << "Size of diffusionCoeffs list not equal to the number of liquid "
            << "species" << nl << exit(FatalError);
    }

    if (liquids_.size() != pSat_.size())
    {
        FatalErrorIn
        (
            "Foam::liquidEvaporation<CloudType>::liquidEvaporation\n"
            "(\n"
            "    const dictionary& dict,\n"
            "    CloudType& cloud\n"
            ")\n"
        )   << "Size of saturation pressure lists not equal to the number of "
            << "liquid species" << nl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::liquidEvaporation<CloudType>::~liquidEvaporation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::liquidEvaporation<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::liquidEvaporation<CloudType>::calculate
(
    const scalar T,
    const scalar d,
    const scalarField& Xc,
    const scalarList& dMassMT,
    const vector& Ur
    const scalar Tc,
    const scalar pc,
    const scalar nuc
    const scalar dt,
) const
{
    scalar dMassTot = 0.0;

    if (T < Tvap_)
    {
        // not reached temperature threshold
        return 0.0;
    }
    else
    {
        // droplet area
        scalar A = mathematicalConstant::pi*sqr(d);

        // universal gas constant
        const scalar R = specie::RR.value();

        // Reynolds number
        scalar Re = mag(Ur)*d/(nuc + ROOTVSMALL);

        // calculate mass transfer of each specie
        forAll(dMassMT, i)
        {
            // Schmidt number
            scalar Sc = nuc/(Dm_[i] + ROOTVSMALL);

            // Sherwood number
            scalar Sh = this->Sh(Re, Sc);

            // mass transfer coefficient [m/s]
            scalar kc = Sh*Dm_[i]/(d + ROOTVSMALL);

            // vapour concentration at droplet surface [kgmol/m3]
            scalar Cs = pSat(i, T)/(R*T);

            // vapour concentration in bulk gas [kgmol/m3]
            scalar Cinf = Xc[i]*pc/(R*Tc);

            // molar flux of vapour [kgmol/m2/s]
            scalar Ni = max(kc*(Cs - Cinf), 0.0);

            // mass transfer
            scalar dm = Ni*A*liquids_.properies()[i].W()*dt;
            dMassMT[i] -= dm;
            dMassTot += dm;
        }
    }

    return dMassTot;
}


// ************************************************************************* //

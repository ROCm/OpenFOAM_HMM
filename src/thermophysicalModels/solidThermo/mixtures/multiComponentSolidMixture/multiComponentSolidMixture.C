/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "multiComponentSolidMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
const ThermoType&  Foam::multiComponentSolidMixture<ThermoType>::
constructSolidData
(
    const dictionary& thermoSolidDict
)
{
    forAll(components_, i)
    {
        solidData_.set
        (
            i,
            new ThermoType
            (
                thermoSolidDict.subDict(components_[i] + "Coeffs")
            )
        );
    }

    return solidData_[0];
}


template<class ThermoType>
void Foam::multiComponentSolidMixture<ThermoType>::correctMassFractions()
{
    volScalarField Yt("Yt", Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }


}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::X
(
    label specieI, label celli, scalar p, scalar T
) const
{
    scalar rhoInv = 0.0;
    forAll(solidData_, i)
    {
        rhoInv += Y_[i][celli]/solidData_[i].rho(p, T);
    }

    scalar X = Y_[specieI][celli]/solidData_[specieI].rho(p, T);

    return (X/rhoInv);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multiComponentSolidMixture<ThermoType>::multiComponentSolidMixture
(
    const dictionary& thermoSolidDict,
    const fvMesh& mesh
)
:
    basicSolidMixture
    (
        thermoSolidDict.lookup("solidComponents"),
        mesh
    ),
    solidData_(components_.size()),
    mixture_("mixture", constructSolidData(thermoSolidDict)),
    mixtureVol_("mixture", constructSolidData(thermoSolidDict))
{
    correctMassFractions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class ThermoType>
const ThermoType& Foam::multiComponentSolidMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(solidData_, i)
    {
        rhoInv += Y_[i][celli]/solidData_[i].rho(p, T);
    }

    mixtureVol_ = Y_[0][celli]/solidData_[0].rho(p, T)/rhoInv*solidData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n][celli]/solidData_[n].rho(p, T)/rhoInv*solidData_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentSolidMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = Y_[0][celli]*solidData_[0];
    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n][celli]*solidData_[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentSolidMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(solidData_, i)
    {
        rhoInv += Y_[i].boundaryField()[patchi][facei]/solidData_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0].boundaryField()[patchi][facei]
      / solidData_[0].rho(p, T)
      / rhoInv
      * solidData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n].boundaryField()[patchi][facei]
          / solidData_[n].rho(p,T)
          / rhoInv
          * solidData_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentSolidMixture<ThermoType>::
patchFaceMixture
(
    const label patchi,
    const label facei
) const
{

    mixture_ =
        Y_[0].boundaryField()[patchi][facei]*solidData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n].boundaryField()[patchi][facei]*solidData_[n];
    }

    return mixture_;
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::Cp
(
    const label specieI, const scalar p, const scalar T
) const
{
    return solidData_[specieI].Cp(p, T);
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::Cv
(
    const label specieI, const scalar p, const scalar T
) const
{
    return solidData_[specieI].Cv(p, T);
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::Ha
(
    const label specieI, const scalar p, const scalar T
) const
{
    return solidData_[specieI].Ha(p, T);
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::Hs
(
    const label specieI, const scalar p, const scalar T
) const
{
    return solidData_[specieI].Hs(p, T);
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::Hc
(
    const label specieI
) const
{
    return solidData_[specieI].Hc();
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::rho
(
    const label specieI, const scalar p, const scalar T
) const
{
    return solidData_[specieI].rho(p, T);
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::kappa
(
    const label specieI, const scalar p, const scalar T
) const
{
    return solidData_[specieI].kappa(T);
}


template<class ThermoType>
Foam::vector Foam::multiComponentSolidMixture<ThermoType>::Kappa
(
    const label specieI, const scalar p, const scalar T
) const
{
    return solidData_[specieI].Kappa(T);
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::alpha
(
    const label specieI, const scalar p, const scalar T
) const
{
    return solidData_[specieI].kappa(T)/solidData_[specieI].Cpv(p, T);
}



template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::rho
(
    scalar p, scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].rho(p, T)*X(i, celli, p, T);
    }
    return tmp;
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::kappaRad
(
    scalar p, scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].kappaRad(T)*X(i, celli, p, T);
    }
    return tmp;
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::sigmaS
(
    scalar p, scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].sigmaS(T)*X(i, celli, p, T);
    }
    return tmp;
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::kappa
(
    scalar p, scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].kappa(T)*X(i, celli, p, T);
    }
    return tmp;
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::emissivity
(
    scalar p, scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].emissivity(T)*X(i, celli, p, T);
    }
    return tmp;
}


template<class ThermoType>
Foam::scalar Foam::multiComponentSolidMixture<ThermoType>::Cp
(
    scalar p, scalar T, label celli
) const
{
    scalar tmp = 0.0;
    forAll(solidData_, i)
    {
        tmp += solidData_[i].Cp(p, T)*Y_[i][celli];
    }
    return tmp;
}


template<class ThermoType>
void Foam::multiComponentSolidMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    forAll(components_, i)
    {
        solidData_[i] =
            ThermoType(thermoDict.subDict(components_[i] + "Coeffs"));
    }
}


// ************************************************************************* //

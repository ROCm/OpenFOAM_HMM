/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "hhuMixtureThermo.H"
#include "fvMesh.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
hhuMixtureThermo<MixtureType>::hhuMixtureThermo(const fvMesh& mesh)
:
    hhuCombustionThermo(mesh),
    MixtureType(*this, mesh)
{
    forAll(h_, celli)
    {
        h_[celli] = this->cellMixture(celli).H(T_[celli]);
        hu_[celli] = this->cellReactants(celli).H(Tu_[celli]);
    }

    forAll(h_.boundaryField(), patchi)
    {
        h_.boundaryField()[patchi] == h(T_.boundaryField()[patchi], patchi);

        fvPatchScalarField& phu = hu_.boundaryField()[patchi];
        const fvPatchScalarField& pTu = Tu_.boundaryField()[patchi];

        forAll(phu, facei)
        {
            phu[facei] = this->patchFaceReactants(patchi, facei).H(pTu[facei]);
        }
    }

    hBoundaryCorrection(h_);
    huBoundaryCorrection(hu_);

    calculate();
    psi_.oldTime();   // Switch on saving old time
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType>
hhuMixtureThermo<MixtureType>::~hhuMixtureThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
void hhuMixtureThermo<MixtureType>::calculate()
{
    forAll(T_, celli)
    {
        const typename MixtureType::thermoType& mixture_ = 
            this->cellMixture(celli);

        T_[celli] = mixture_.TH(h_[celli], T_[celli]);
        psi_[celli] = mixture_.psi(p_[celli], T_[celli]);

        mu_[celli] = mixture_.mu(T_[celli]);
        alpha_[celli] = mixture_.alpha(T_[celli]);

        Tu_[celli] = this->cellReactants(celli).TH(hu_[celli], Tu_[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = p_.boundaryField()[patchi];
        fvPatchScalarField& pT = T_.boundaryField()[patchi];
        fvPatchScalarField& pTu = Tu_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = psi_.boundaryField()[patchi];

        fvPatchScalarField& ph = h_.boundaryField()[patchi];
        fvPatchScalarField& phu = hu_.boundaryField()[patchi];

        fvPatchScalarField& pmu_ = mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha_ = alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = mixture_.H(pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu_[facei] = mixture_.mu(pT[facei]);
                palpha_[facei] = mixture_.alpha(pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.TH(ph[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu_[facei] = mixture_.mu(pT[facei]);
                palpha_[facei] = mixture_.alpha(pT[facei]);

                pTu[facei] =
                    this->patchFaceReactants(patchi, facei)
                   .TH(phu[facei], pTu[facei]);
            }
        }
    }
}


template<class MixtureType>
void hhuMixtureThermo<MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering hhuMixtureThermo<MixtureType>::correct()" << endl;
    }

    // force the saving of the old-time values
    psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "exiting hhuMixtureThermo<MixtureType>::correct()" << endl;
    }
}

template<class MixtureType>
tmp<scalarField> hhuMixtureThermo<MixtureType>::h
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    forAll(T, celli)
    {
        h[celli] = this->cellMixture(cells[celli]).H(T[celli]);
    }

    return th;
}


template<class MixtureType>
tmp<scalarField> hhuMixtureThermo<MixtureType>::h
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    forAll(T, facei)
    {
        h[facei] = this->patchFaceMixture(patchi, facei).H(T[facei]);
    }

    return th;
}


template<class MixtureType>
tmp<scalarField> hhuMixtureThermo<MixtureType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));

    scalarField& cp = tCp();

    forAll(T, facei)
    {
        cp[facei] = this->patchFaceMixture(patchi, facei).Cp(T[facei]);
    }

    return tCp;
}


template<class MixtureType>
tmp<volScalarField> hhuMixtureThermo<MixtureType>::Cp() const
{
    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );

    volScalarField& cp = tCp();

    forAll(T_, celli)
    {
        cp[celli] = this->cellMixture(celli).Cp(T_[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        cp.boundaryField()[patchi] = Cp(T_.boundaryField()[patchi], patchi);
    }

    return tCp;
}



template<class MixtureType>
tmp<scalarField> hhuMixtureThermo<MixtureType>::hu
(
    const scalarField& Tu,
    const labelList& cells
) const
{
    tmp<scalarField> thu(new scalarField(Tu.size()));
    scalarField& hu = thu();

    forAll(Tu, celli)
    {
        hu[celli] = this->cellReactants(cells[celli]).H(Tu[celli]);
    }

    return thu;
}


template<class MixtureType>
tmp<scalarField> hhuMixtureThermo<MixtureType>::hu
(
    const scalarField& Tu,
    const label patchi
) const
{
    tmp<scalarField> thu(new scalarField(Tu.size()));
    scalarField& hu = thu();

    forAll(Tu, facei)
    {
        hu[facei] = this->patchFaceReactants(patchi, facei).H(Tu[facei]);
    }

    return thu;
}


template<class MixtureType>
tmp<volScalarField> hhuMixtureThermo<MixtureType>::Tb() const
{
    tmp<volScalarField> tTb
    (
        new volScalarField
        (
            IOobject
            (
                "Tb",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_
        )
    );

    volScalarField& Tb_ = tTb();

    forAll(Tb_, celli)
    {
        Tb_[celli] = this->cellProducts(celli).TH(h_[celli], T_[celli]);
    }

    forAll(Tb_.boundaryField(), patchi)
    {
        fvPatchScalarField& pTb = Tb_.boundaryField()[patchi];

        const fvPatchScalarField& ph = h_.boundaryField()[patchi];
        const fvPatchScalarField& pT = T_.boundaryField()[patchi];

        forAll(pTb, facei)
        {
            pTb[facei] =
                this->patchFaceProducts(patchi, facei).TH(ph[facei], pT[facei]);
        }
    }

    return tTb;
}


template<class MixtureType>
tmp<volScalarField> hhuMixtureThermo<MixtureType>::psiu() const
{
    tmp<volScalarField> tpsiu
    (
        new volScalarField
        (
            IOobject
            (
                "psiu",
                psi_.time().timeName(),
                psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            psi_.dimensions()
        )
    );

    volScalarField& psiu = tpsiu();

    forAll(psiu, celli)
    {
        psiu[celli] = this->cellReactants(celli).psi(p_[celli], Tu_[celli]);
    }

    forAll(psiu.boundaryField(), patchi)
    {
        fvPatchScalarField& ppsiu = psiu.boundaryField()[patchi];

        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const fvPatchScalarField& pTu = Tu_.boundaryField()[patchi];

        forAll(ppsiu, facei)
        {
            ppsiu[facei] =
                this->
                patchFaceReactants(patchi, facei).psi(pp[facei], pTu[facei]);
        }
    }

    return tpsiu;
}


template<class MixtureType>
tmp<volScalarField> hhuMixtureThermo<MixtureType>::psib() const
{
    tmp<volScalarField> tpsib
    (
        new volScalarField
        (
            IOobject
            (
                "psib",
                psi_.time().timeName(),
                psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            psi_.dimensions()
        )
    );

    volScalarField& psib = tpsib();

    volScalarField Tb_ = Tb();

    forAll(psib, celli)
    {
        psib[celli] = this->cellReactants(celli).psi(p_[celli], Tb_[celli]);
    }

    forAll(psib.boundaryField(), patchi)
    {
        fvPatchScalarField& ppsib = psib.boundaryField()[patchi];

        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const fvPatchScalarField& pTb = Tb_.boundaryField()[patchi];

        forAll(ppsib, facei)
        {
            ppsib[facei] =
                this->patchFaceReactants(patchi, facei).psi(pp[facei], pTb[facei]);
        }
    }

    return tpsib;
}


template<class MixtureType>
tmp<volScalarField> hhuMixtureThermo<MixtureType>::muu() const
{
    tmp<volScalarField> tmuu
    (
        new volScalarField
        (
            IOobject
            (
                "muu",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimensionSet(1, -1, -1, 0, 0)
        )
    );

    volScalarField& muu_ = tmuu();

    forAll(muu_, celli)
    {
        muu_[celli] = this->cellReactants(celli).mu(Tu_[celli]);
    }

    forAll(muu_.boundaryField(), patchi)
    {
        fvPatchScalarField& pMuu = muu_.boundaryField()[patchi];
        const fvPatchScalarField& pTu = Tu_.boundaryField()[patchi];

        forAll(pMuu, facei)
        {
            pMuu[facei] =
                this->patchFaceReactants(patchi, facei).mu(pTu[facei]);
        }
    }

    return tmuu;
}


template<class MixtureType>
tmp<volScalarField> hhuMixtureThermo<MixtureType>::mub() const
{
    tmp<volScalarField> tmub
    (
        new volScalarField
        (
            IOobject
            (
                "mub",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimensionSet(1, -1, -1, 0, 0)
        )
    );

    volScalarField& mub_ = tmub();
    volScalarField Tb_ = Tb();

    forAll(mub_, celli)
    {
        mub_[celli] = this->cellProducts(celli).mu(Tb_[celli]);
    }

    forAll(mub_.boundaryField(), patchi)
    {
        fvPatchScalarField& pMub = mub_.boundaryField()[patchi];
        const fvPatchScalarField& pTb = Tb_.boundaryField()[patchi];

        forAll(pMub, facei)
        {
            pMub[facei] = this->patchFaceProducts(patchi, facei).mu(pTb[facei]);
        }
    }

    return tmub;
}


template<class MixtureType>
bool hhuMixtureThermo<MixtureType>::read()
{
    if (hhuCombustionThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

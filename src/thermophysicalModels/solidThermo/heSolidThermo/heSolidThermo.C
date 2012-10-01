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

#include "heSolidThermo.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
void Foam::heSolidThermo<BasicSolidThermo, MixtureType>::calculate()
{
    scalarField& TCells = this->T_.internalField();

    const scalarField& hCells = this->he_.internalField();
    const scalarField& pCells = this->p_.internalField();
    scalarField& rhoCells = this->rho_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        const typename MixtureType::thermoType& volMixture_ =
            this->cellVolMixture(pCells[celli], TCells[celli], celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        rhoCells[celli] = volMixture_.rho(pCells[celli], TCells[celli]);

        alphaCells[celli] =
            volMixture_.kappa(TCells[celli])
            /
            mixture_.Cpv(pCells[celli], TCells[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];

        fvPatchScalarField& ph = this->he_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                const typename MixtureType::thermoType& volMixture_ =
                    this->patchFaceVolMixture
                    (
                        pp[facei],
                        pT[facei],
                        patchi,
                        facei
                    );

                ph[facei] = mixture_.HE(pp[facei], pT[facei]);
                prho[facei] = volMixture_.rho(pp[facei], pT[facei]);
                palpha[facei] =
                    volMixture_.kappa(pT[facei])
                  / mixture_.Cpv(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                const typename MixtureType::thermoType& volMixture_ =
                    this->patchFaceVolMixture
                    (
                        pp[facei],
                        pT[facei],
                        patchi,
                        facei
                    );

                pT[facei] = mixture_.THE(ph[facei], pp[facei] ,pT[facei]);
                prho[facei] = volMixture_.rho(pp[facei], pT[facei]);
                palpha[facei] =
                    volMixture_.kappa(pT[facei])
                  / mixture_.Cpv(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::
heSolidThermo(const fvMesh& mesh)
:
    heThermo<BasicSolidThermo, MixtureType>(mesh)
{
    calculate();
}


template<class BasicSolidThermo, class MixtureType>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::
heSolidThermo(const fvMesh& mesh, const dictionary& dict)
:
    heThermo<BasicSolidThermo, MixtureType>(mesh, dict)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::~heSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
void Foam::heSolidThermo<BasicSolidThermo, MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering heSolidThermo<MixtureType>::correct()" << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "exiting heSolidThermo<MixtureType>::correct()" << endl;
    }
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::volVectorField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::Kappa() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volVectorField> tKappa
    (
        new volVectorField
        (
            IOobject
            (
                "Kappa",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimTime/dimLength/dimTemperature
        )
    );

    volVectorField& Kappa = tKappa();
    vectorField& KappaCells = Kappa.internalField();
    const scalarField& TCells = this->T_.internalField();
    const scalarField& pCells = this->p_.internalField();

    forAll(KappaCells, celli)
    {
        Kappa[celli] =
            this->cellVolMixture
            (
                pCells[celli],
                TCells[celli],
                celli
            ).Kappa(TCells[celli]);
    }

    forAll(Kappa.boundaryField(), patchi)
    {
        vectorField& Kappap = Kappa.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        const scalarField& pp = this->p_.boundaryField()[patchi];

        forAll(Kappap, facei)
        {
            Kappap[facei] =
                this->patchFaceVolMixture
                (
                    pp[facei],
                    pT[facei],
                    patchi,
                    facei
                ).Kappa(pT[facei]);
        }
    }

    return tKappa;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::kappaRad() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tkappaRad
    (
        new volScalarField
        (
            IOobject
            (
                "kappaRad",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            inv(dimLength)
        )
    );

    volScalarField& kappaRad = tkappaRad();
    scalarField& kappaRadCells = kappaRad.internalField();
    const scalarField& TCells = this->T_.internalField();
    const scalarField& pCells = this->p_.internalField();

    forAll(kappaRadCells, celli)
    {
        kappaRadCells[celli] =
            this->cellVolMixture
            (
                pCells[celli],
                TCells[celli],
                celli
            ).kappaRad(TCells[celli]);
    }

    forAll(kappaRad.boundaryField(), patchi)
    {
        scalarField& kappaRadp = kappaRad.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        const scalarField& pp = this->p_.boundaryField()[patchi];

        forAll(kappaRadp, facei)
        {
            kappaRadp[facei] =
                this->patchFaceVolMixture
                (
                    pp[facei],
                    pT[facei],
                    patchi,
                    facei
                ).kappaRad(pT[facei]);
        }
    }

    return tkappaRad;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::sigmaS() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tsigmaS
    (
        new volScalarField
        (
            IOobject
            (
                "sigmaS",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            inv(dimLength)
        )
    );

    volScalarField& sigmaS = tsigmaS();
    scalarField& sigmaSCells = sigmaS.internalField();
    const scalarField& TCells = this->T_.internalField();
    const scalarField& pCells = this->p_.internalField();

    forAll(sigmaSCells, celli)
    {
        sigmaSCells[celli] =
            this->cellVolMixture
            (
                pCells[celli],
                TCells[celli],
                celli
            ).sigmaS(TCells[celli]);
    }

    forAll(sigmaS.boundaryField(), patchi)
    {
        scalarField& sigmaSp = sigmaS.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        const scalarField& pp = this->p_.boundaryField()[patchi];

        forAll(sigmaSp, facei)
        {
            sigmaSp[facei] =
                this->patchFaceVolMixture
                (
                    pp[facei],
                    pT[facei],
                    patchi,
                    facei
                ).sigmaS(pT[facei]);
        }
    }

    return tsigmaS;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::emissivity() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> temissivity
    (
        new volScalarField
        (
            IOobject
            (
                "emissivity",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            inv(dimLength)
        )
    );

    volScalarField& emissivity = temissivity();
    scalarField& emissivityCells = emissivity.internalField();
    const scalarField& TCells = this->T_.internalField();
    const scalarField& pCells = this->p_.internalField();

    forAll(emissivityCells, celli)
    {
        emissivityCells[celli] =
            this->cellVolMixture
            (
                pCells[celli],
                TCells[celli],
                celli
            ).emissivity(TCells[celli]);
    }

    forAll(emissivity.boundaryField(), patchi)
    {
        scalarField& emissivityp = emissivity.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        const scalarField& pp = this->p_.boundaryField()[patchi];

        forAll(emissivityp, facei)
        {
            emissivityp[facei] =
                this->patchFaceVolMixture
                (
                    pp[facei],
                    pT[facei],
                    patchi,
                    facei
                ).emissivity(pT[facei]);
        }
    }

    return temissivity;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::vectorField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::Kappa
(
    const label patchi
) const
{
    const scalarField& pp = this->p_.boundaryField()[patchi];
    const scalarField& Tp = this->T_.boundaryField()[patchi];
    tmp<vectorField> tKappa(new vectorField(pp.size()));

    vectorField& Kappap = tKappa();

    forAll(Tp, facei)
    {
        Kappap[facei] =
            this->patchFaceVolMixture
            (
                pp[facei],
                Tp[facei],
                patchi,
                facei
            ).Kappa(Tp[facei]);
    }

    return tKappa;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::kappaRad
(
    const label patchi
) const
{
    const scalarField& Tp = this->T_.boundaryField()[patchi];
    tmp<scalarField> tKappaRad(new scalarField(Tp.size()));
    scalarField& KappapRadp = tKappaRad();
    const scalarField& pp = this->p_.boundaryField()[patchi];

    forAll(Tp, facei)
    {
        KappapRadp[facei] =
            this->patchFaceVolMixture
            (
                pp[facei],
                Tp[facei],
                patchi,
                facei
            ).kappaRad(Tp[facei]);
    }

    return tKappaRad;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::sigmaS
(
    const label patchi
) const
{
    const scalarField& Tp = this->T_.boundaryField()[patchi];
    tmp<scalarField> tsigmaS(new scalarField(Tp.size()));
    scalarField& sigmaSp = tsigmaS();
    const scalarField& pp = this->p_.boundaryField()[patchi];


    forAll(Tp, facei)
    {
        sigmaSp[facei] =
            this->patchFaceVolMixture
            (
                pp[facei],
                Tp[facei],
                patchi,
                facei
            ).sigmaS(Tp[facei]);
    }

    return tsigmaS;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::emissivity
(
    const label patchi
) const
{
    const scalarField& Tp = this->T_.boundaryField()[patchi];
    tmp<scalarField> temissivity(new scalarField(Tp.size()));
    scalarField& emissivity = temissivity();
    const scalarField& pp = this->p_.boundaryField()[patchi];

    forAll(Tp, facei)
    {
        emissivity[facei] =
            this->patchFaceVolMixture
            (
                pp[facei],
                Tp[facei],
                patchi,
                facei
            ).emissivity(Tp[facei]);
    }

    return temissivity;
}


// ************************************************************************* //

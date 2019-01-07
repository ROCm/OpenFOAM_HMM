/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "MassTransferPhaseSystem.H"

#include "HashPtrTable.H"

#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvMatrix.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::MassTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels("massTransferModel", massTransferModels_);

    forAllConstIters(massTransferModels_, iterModel)
    {
        if (!dmdt_.found(iterModel()->pair()))
        {
            dmdt_.set
            (
                iterModel()->pair(),
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("dmdt",iterModel()->pair().name()),
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, Zero)
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::calculateL
(
    const volScalarField& dmdtNetki,
    const phasePairKey& keyik,
    const phasePairKey& keyki,
    const volScalarField& T
) const
{
    tmp<volScalarField> tL
    (
        new volScalarField
        (
            IOobject
            (
                "tL",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimMass, Zero)
        )
    );
    volScalarField& L = tL.ref();

    if (massTransferModels_.found(keyik))
    {
        const autoPtr<interfaceCompositionModel>& interfacePtr =
            massTransferModels_[keyik];

        word speciesName = interfacePtr->transferSpecie();

        auto tempOpen = speciesName.find('.');

        //const word species(speciesName(0, tempOpen));
        const word species(speciesName.substr(0, tempOpen));

        L -= neg(dmdtNetki)*interfacePtr->L(species, T);
    }

    if (massTransferModels_.found(keyki))
    {
        const autoPtr<interfaceCompositionModel>& interfacePtr =
            massTransferModels_[keyki];

        word speciesName = interfacePtr->transferSpecie();

        auto tempOpen = speciesName.find('.');

        const word species(speciesName.substr(0, tempOpen));

        L += pos(dmdtNetki)*interfacePtr->L(species, T);
    }

    return tL;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    tmp<volScalarField> tdmdt
    (
        new volScalarField
        (
            IOobject
            (
                "dmdt",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimDensity/dimTime, Zero)
        )
    );

    volScalarField& dmdt = tdmdt.ref();

    if (dmdt_.found(key))
    {
        dmdt = *dmdt_[key];
    }

    return tdmdt;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::heatTransfer
(
    const volScalarField& T
)
{
    tmp<fvScalarMatrix> tEqnPtr
    (
        new fvScalarMatrix(T, dimEnergy/dimTime)
    );

    fvScalarMatrix& eqn = tEqnPtr.ref();

    forAllConstIters(this->phaseModels_, iteri)
    {
        const phaseModel& phasei = iteri()();

        auto iterk = iteri;

        for (++iterk; iterk != this->phaseModels_.end(); ++iterk)
        {
            if (iteri()().name() != iterk()().name())
            {
                const phaseModel& phasek = iterk()();

                // Phase i to phase k
                const phasePairKey keyik(phasei.name(), phasek.name(), true);

                // Phase k to phase i
                const phasePairKey keyki(phasek.name(), phasei.name(), true);

                // Net mass transfer from k to i phase
                tmp<volScalarField> tdmdtNetki
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "tdmdtYki",
                            this->mesh().time().timeName(),
                            this->mesh(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        this->mesh(),
                        dimensionedScalar(dimDensity/dimTime, Zero)
                    )
                );
                volScalarField& dmdtNetki = tdmdtNetki.ref();


                if (massTransferModels_.found(keyik))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyik];

                    // Explicit temperature mass transfer rate
                    tmp<volScalarField> Kexp =
                        interfacePtr->Kexp(interfaceCompositionModel::T, T);

                    if (Kexp.valid())
                    {
                        dmdtNetki -= Kexp.ref();
                        *dmdt_[keyik] = Kexp.ref();
                    }
                }

                // Looking for mass transfer in the other direction (k to i)
                if (massTransferModels_.found(keyki))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyki];

                    // Explicit temperature mass transfer rate
                    const tmp<volScalarField> Kexp =
                        interfacePtr->Kexp(interfaceCompositionModel::T, T);

                    if (Kexp.valid())
                    {
                        dmdtNetki += Kexp.ref();
                        *dmdt_[keyki] = Kexp.ref();
                    }

                }

                word keyikName(phasei.name() + phasek.name());
                word keykiName(phasek.name() + phasei.name());

                eqn -=
                    (
                        dmdtNetki
                        *(
                            calculateL(dmdtNetki, keyik, keyki, T)
                            - (phasek.Cp() - phasei.Cp())
                            * constant::standard::Tstd
                        )
                    );
            }
        }
    }
    return tEqnPtr;
}


template<class BasePhaseSystem>
void Foam::MassTransferPhaseSystem<BasePhaseSystem>::massSpeciesTransfer
(
    const phaseModel& phase,
    volScalarField::Internal& Su,
    volScalarField::Internal& Sp,
    const word speciesName
)
{
    // Fill the volumetric mass transfer for species
    forAllConstIters(massTransferModels_, iter)
    {
        if (iter()->transferSpecie() == speciesName)
        {
            // Explicit source
            Su +=
                  this->Su()[phase.name()]
                + this->Sp()[phase.name()]*phase.oldTime();
        }
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2015-2017 OpenFOAM Foundation
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

#include "phaseSystem.H"
#include "surfaceTensionModel.H"
#include "aspectRatioModel.H"
#include "surfaceInterpolate.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
}

const Foam::word Foam::phaseSystem::propertiesName("phaseProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::calcPhi
(
    const phaseModelList& phaseModels
) const
{
    tmp<surfaceScalarField> tmpPhi
    (
        new surfaceScalarField
        (
            "phi",
            fvc::interpolate(phaseModels[0])*phaseModels[0].phi()
        )
    );

    for (label phasei=1; phasei<phaseModels.size(); ++phasei)
    {
        tmpPhi.ref() +=
            fvc::interpolate(phaseModels[phasei])*phaseModels[phasei].phi();
    }

    return tmpPhi;
}


void Foam::phaseSystem::generatePairs
(
    const dictTable& modelDicts
)
{
    forAllConstIters(modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        // pair already exists
        if (phasePairs_.found(key))
        {}

        // new ordered pair
        else if (key.ordered())
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new orderedPhasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }

        // new unordered pair
        else
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new phasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    phaseModels_(lookup("phases"), phaseModel::iNew(*this)),

    phi_(calcPhi(phaseModels_)),

    dpdt_
    (
        IOobject
        (
            "dpdt",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimPressure/dimTime, Zero)
    ),

    MRF_(mesh_)
{
    phi_.writeOpt() = IOobject::AUTO_WRITE;

    // Blending methods
    for (const entry& dEntry : subDict("blending"))
    {
        blendingMethods_.insert
        (
            dEntry.dict().dictName(),
            blendingMethod::New
            (
                dEntry.dict(),
                phaseModels_.toc()
            )
        );
    }

    // Sub-models
    generatePairsAndSubModels("surfaceTension", surfaceTensionModels_);
    generatePairsAndSubModels("aspectRatio", aspectRatioModels_);

    correctKinematics();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseSystem::rho() const
{
    auto phasei = phaseModels_.cbegin();

    tmp<volScalarField> tmpRho
    (
        phasei() * phasei().rho()
    );

    for (++phasei; phasei != phaseModels_.cend(); ++phasei)
    {
        tmpRho.ref() += phasei() * phasei().rho();
    }

    return tmpRho;
}


Foam::tmp<Foam::volVectorField> Foam::phaseSystem::U() const
{
    auto phasei = phaseModels_.cbegin();

    tmp<volVectorField> tmpU
    (
        phasei() * phasei().U()
    );

    for (++phasei; phasei != phaseModels_.cend(); ++phasei)
    {
        tmpU.ref() += phasei() * phasei().U();
    }

    return tmpU;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::E(const phasePairKey& key) const
{
    if (aspectRatioModels_.found(key))
    {
        return aspectRatioModels_[key]->E();
    }

    return tmp<volScalarField>::New
    (
        IOobject
        (
            aspectRatioModel::typeName + ":E",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar("one", dimless, 1)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::sigma(const phasePairKey& key) const
{
    if (surfaceTensionModels_.found(key))
    {
        return surfaceTensionModels_[key]->sigma();
    }

    return tmp<volScalarField>::New
    (
        IOobject
        (
            surfaceTensionModel::typeName + ":sigma",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(surfaceTensionModel::dimSigma, Zero)
    );
}


void Foam::phaseSystem::solve()
{}


void Foam::phaseSystem::correct()
{
    for (phaseModel& phase : phaseModels_)
    {
        phase.correct();
    }
}


void Foam::phaseSystem::correctKinematics()
{
    bool updateDpdt = false;

    for (phaseModel& phase : phaseModels_)
    {
        phase.correctKinematics();

        updateDpdt = updateDpdt || phase.thermo().dpdt();
    }

    // Update the pressure time-derivative if required
    if (updateDpdt)
    {
        dpdt_ = fvc::ddt(phaseModels_.cbegin()().thermo().p());
    }
}


void Foam::phaseSystem::correctThermo()
{
    for (phaseModel& phase : phaseModels_)
    {
        phase.correctThermo();
    }
}


void Foam::phaseSystem::correctTurbulence()
{
    for (phaseModel& phase : phaseModels_)
    {
        phase.correctTurbulence();
    }
}


void Foam::phaseSystem::correctEnergyTransport()
{
    for (phaseModel& phase : phaseModels_)
    {
        phase.correctEnergyTransport();
    }
}


bool Foam::phaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        for (phaseModel& phase : phaseModels_)
        {
            readOK &= phase.read();
        }

        // models ...

        return readOK;
    }

    return false;
}


// ************************************************************************* //

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

#include "phaseSystem.H"
#include "surfaceTensionModel.H"
#include "porousModel.H"

#include "HashPtrTable.H"

#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvMatrix.H"

#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
}

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

const Foam::word Foam::phaseSystem::phasePropertiesName("phaseProperties");

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::phaseSystem::phaseModelTable
Foam::phaseSystem::generatePhaseModels(const wordList& phaseNames) const
{
    phaseModelTable phaseModels;

    forAllConstIter(wordList, phaseNames, phaseNameIter)
    {
        phaseModels.insert
        (
            *phaseNameIter,
            phaseModel::New
            (
                *this,
                *phaseNameIter
            )
        );
    }

    return phaseModels;
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::generatePhi
(
    const phaseModelTable& phaseModels
) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels.begin();

    auto tmpPhi = tmp<surfaceScalarField>::New
    (
        "phi",
        fvc::interpolate(phaseModelIter()())*phaseModelIter()->phi()
    );

    ++phaseModelIter;

    for (; phaseModelIter != phaseModels.end(); ++phaseModelIter)
    {
        tmpPhi.ref() +=
            fvc::interpolate(phaseModelIter()())
           *phaseModelIter()->phi();
    }

    return tmpPhi;
}


void Foam::phaseSystem::generatePairs(const dictTable& modelDicts)
{
    forAllConstIter(dictTable, modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        // pair already exists
        if (phasePairs_.found(key))
        {
            // do nothing ...
        }
        else if (key.ordered())
        {
            // New ordered pair
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
        else
        {
            // New unordered pair
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


void Foam::phaseSystem::generatePairsTable()
{

    forAllConstIter(phaseModelTable, phaseModels_, phaseIter1)
    {
        forAllConstIter(phaseModelTable, phaseModels_, phaseIter2)
        {
            if (phaseIter1()->name() != phaseIter2()->name())
            {
                phasePairKey key
                (
                    phaseIter1()->name(),
                    phaseIter2()->name(),
                    true
                );

                phasePairKey keyInverse
                (
                    phaseIter2()->name(),
                    phaseIter1()->name(),
                    true
                );

                if
                (
                    !totalPhasePairs_.found(key)
                 && !totalPhasePairs_.found(keyInverse)
                )
                {
                    totalPhasePairs_.set
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
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    basicThermo(mesh, word::null, phasePropertiesName),
    mesh_(mesh),
    phaseNames_(lookup("phases")),
    phi_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVolume/dimTime, Zero)
    ),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimTime, Zero)
    ),
    phaseModels_(generatePhaseModels(phaseNames_)),
    phasePairs_(),
    totalPhasePairs_(),
    Prt_
    (   dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt", *this, 1.0
        )
    )
{
    rhoPhi_.setOriented();
    phi_.setOriented();

    // sub models
    if (found("surfaceTension"))
    {
        generatePairsAndSubModels
        (
            "surfaceTension",
            mesh_,
            surfaceTensionModels_
        );
    }
    if (found("interfacePorous"))
    {
        generatePairsAndSubModels
        (
            "interfacePorous",
            mesh_,
            interfacePorousModelTable_
        );
    }

    // Total phase pair
    generatePairsTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseSystem::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::hc() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tAlphaHc
    (
        phaseModelIter()()*phaseModelIter()->hc()
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); phaseModelIter++)
    {
        tAlphaHc.ref() += phaseModelIter()()*phaseModelIter()->hc();
    }

    return tAlphaHc;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::THE
(
    const scalarField& e,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::THE
(
    const scalarField& e,
    const scalarField& p,
    const scalarField& T0,
    const label patchI
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::rho() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tmpRho
    (
        phaseModelIter()()*phaseModelIter()->rho()
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); phaseModelIter++)
    {
        tmpRho.ref() += phaseModelIter()()*phaseModelIter()->rho();
    }

    return tmpRho;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::rho(const label patchI) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<scalarField> tmpRho
    (
        phaseModelIter()().boundaryField()[patchI]
      * phaseModelIter()->rho()().boundaryField()[patchI]
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); phaseModelIter++)
    {
        tmpRho.ref() +=
            phaseModelIter()().boundaryField()[patchI]
          * phaseModelIter()->rho()().boundaryField()[patchI];
    }

    return tmpRho;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::Cp() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tmpCp
    (
        phaseModelIter()()*phaseModelIter()->Cp()
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); phaseModelIter++)
    {
        tmpCp.ref() += phaseModelIter()()*phaseModelIter()->Cp();
    }

    return tmpCp;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<scalarField> tmpCp
    (
        phaseModelIter()()*phaseModelIter()->Cp(p, T, patchI)
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++phaseModelIter)
    {
        tmpCp.ref() += phaseModelIter()()*phaseModelIter()->Cp(p, T, patchI);
    }

    return tmpCp;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::Cv() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tmpCv
    (
        phaseModelIter()()*phaseModelIter()->Cv()
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); phaseModelIter++)
    {
        tmpCv.ref() += phaseModelIter()()*phaseModelIter()->Cv();
    }

    return tmpCv;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<scalarField> tmpCv
    (
        phaseModelIter()()*phaseModelIter()->Cv(p, T, patchI)
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpCv.ref() += phaseModelIter()()*phaseModelIter()->Cv(p, T, patchI);
    }

    return tmpCv;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::gamma() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tmpCp
    (
        phaseModelIter()()*phaseModelIter()->Cp()
    );

    tmp<volScalarField> tmpCv
    (
        phaseModelIter()()*phaseModelIter()->Cv()
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); phaseModelIter++)
    {
        tmpCp.ref() += phaseModelIter()()*phaseModelIter()->Cp();
        tmpCv.ref() += phaseModelIter()()*phaseModelIter()->Cv();
    }

    return (tmpCp/tmpCv);
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    return
    (
        gamma()().boundaryField()[patchI]
    );
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::Cpv() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tmpCpv
    (
        phaseModelIter()()*phaseModelIter()->Cpv()
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); phaseModelIter++)
    {
        tmpCpv.ref() += phaseModelIter()()*phaseModelIter()->Cpv();
    }

    return tmpCpv;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<scalarField> tmpCpv
    (
        phaseModelIter()()*phaseModelIter()->Cpv(p, T, patchI)
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpCpv.ref() += phaseModelIter()()*phaseModelIter()->Cpv(p, T, patchI);
    }

    return tmpCpv;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::CpByCpv() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tmpCpByCpv
    (
        phaseModelIter()()*phaseModelIter()->CpByCpv()
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); phaseModelIter++)
    {
        tmpCpByCpv.ref() += phaseModelIter()()*phaseModelIter()->CpByCpv();
    }

    return tmpCpByCpv;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchI
) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<scalarField> tmpCpv
    (
        phaseModelIter()().boundaryField()[patchI]
       *phaseModelIter()->CpByCpv(p, T, patchI)
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpCpv.ref() +=
            phaseModelIter()().boundaryField()[patchI]
           *phaseModelIter()->CpByCpv(p, T, patchI);
    }

    return tmpCpv;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::W() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::kappa() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tmpkappa
    (
        phaseModelIter()()*phaseModelIter()->kappa()
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpkappa.ref() += phaseModelIter()()*phaseModelIter()->kappa();
    }

    return tmpkappa;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::kappa(const label patchI) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<scalarField> tmpKappa
    (
        phaseModelIter()().boundaryField()[patchI]
       *phaseModelIter()->kappa(patchI)
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpKappa.ref() +=
            phaseModelIter()().boundaryField()[patchI]
           *phaseModelIter()->kappa(patchI);
    }

    return tmpKappa;
}


Foam::tmp<Foam::volScalarField>Foam::phaseSystem::kappaEff
(
    const volScalarField& kappat
) const
{
    tmp<volScalarField> kappaEff(kappa() + kappat);
    kappaEff.ref().rename("kappaEff");
    return kappaEff;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::kappaEff
(
    const scalarField& kappat,
    const label patchI
) const
{
    return kappa(patchI) + kappat;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::alphaEff
(
    const volScalarField& alphat
) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tmpAlpha
    (
        phaseModelIter()()*phaseModelIter()->alpha()
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpAlpha.ref() += phaseModelIter()()*phaseModelIter()->alpha();
    }

    tmpAlpha.ref() += alphat;

    return tmpAlpha;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::alphaEff
(
    const scalarField& alphat,
    const label patchI
) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<scalarField> tmpAlpha
    (
        phaseModelIter()().boundaryField()[patchI]
       *phaseModelIter()->alpha(patchI)
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpAlpha.ref() +=
            phaseModelIter()().boundaryField()[patchI]
           *phaseModelIter()->alpha(patchI);
    }

    tmpAlpha.ref() += alphat;

    return tmpAlpha;
}


const Foam::dimensionedScalar& Foam::phaseSystem::Prt() const
{
    return Prt_;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::mu() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tmpMu
    (
        phaseModelIter()()*phaseModelIter()->mu()
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpMu.ref() += phaseModelIter()()*phaseModelIter()->mu();
    }

    return tmpMu;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::mu(const label patchI) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<scalarField> tmpMu
    (
        phaseModelIter()().boundaryField()[patchI]
       *phaseModelIter()->mu(patchI)
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpMu.ref() +=
            phaseModelIter()().boundaryField()[patchI]
           *phaseModelIter()->mu(patchI);
    }

    return tmpMu;
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::nu() const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<volScalarField> tmpNu
    (
        phaseModelIter()()*phaseModelIter()->nu()
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpNu.ref() += phaseModelIter()()*phaseModelIter()->nu();
    }

    return tmpNu;
}


Foam::tmp<Foam::scalarField> Foam::phaseSystem::nu(const label patchI) const
{
    phaseModelTable::const_iterator phaseModelIter = phaseModels_.begin();

    tmp<scalarField> tmpNu
    (
        phaseModelIter()().boundaryField()[patchI]
       *phaseModelIter()->nu(patchI)
    );

    ++phaseModelIter;
    for (; phaseModelIter != phaseModels_.end(); ++ phaseModelIter)
    {
        tmpNu.ref() +=
            phaseModelIter()().boundaryField()[patchI]
           *phaseModelIter()->nu(patchI);
    }

    return tmpNu;
}

const Foam::surfaceScalarField& Foam::phaseSystem::phi() const
{
    return phi_;
}


Foam::surfaceScalarField& Foam::phaseSystem::phi()
{
    return phi_;
}


const Foam::surfaceScalarField& Foam::phaseSystem::rhoPhi() const
{
    return rhoPhi_;
}


Foam::surfaceScalarField& Foam::phaseSystem::rhoPhi()
{
    return rhoPhi_;
}


void Foam::phaseSystem::correct()
{
    forAllIter(phaseModelTable, phaseModels_, phaseModelIter)
    {
        phaseModelIter()->correct();
    }
}


void Foam::phaseSystem::correctTurbulence()
{
    forAllIter(phaseModelTable, phaseModels_, phaseModelIter)
    {
        phaseModelIter()->correctTurbulence();
    }
}


const Foam::phaseSystem::phaseModelTable& Foam::phaseSystem::phases() const
{
    return phaseModels_;
}


Foam::phaseSystem::phaseModelTable& Foam::phaseSystem::phases()
{
    return phaseModels_;
}


const Foam::phaseSystem::phasePairTable&
Foam::phaseSystem::totalPhasePairs() const
{
    return totalPhasePairs_;
}


Foam::phaseSystem::phasePairTable& Foam::phaseSystem::totalPhasePairs()
{
    return totalPhasePairs_;
}


bool Foam::phaseSystem::incompressible() const
{
    bool incompressible = true;

    forAllConstIter(phaseModelTable, phaseModels_, phaseModelIter)
    {
        incompressible *= phaseModelIter()->thermo().incompressible();
    }

    return incompressible;
}


bool Foam::phaseSystem::incompressible(const word phaseName) const
{
    return phaseModels_[phaseName]->thermo().incompressible();
}


bool Foam::phaseSystem::isochoric() const
{
     bool isochoric = true;

    forAllConstIter(phaseModelTable, phaseModels_, phaseModelIter)
    {
        isochoric *= phaseModelIter()->thermo().isochoric();
    }

    return isochoric;
}


const Foam::fvMesh& Foam::phaseSystem::mesh() const
{
    return mesh_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::phaseSystem::surfaceTensionForce() const
{
    auto tstf = tmp<surfaceScalarField>::New
    (
        IOobject
        (
            "surfaceTensionForce",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), Zero)
    );

    auto& stf = tstf.ref();
    stf.setOriented();

    if (surfaceTensionModels_.size() > 0)
    {
        forAllConstIter(phaseModelTable, phaseModels_, iter1)
        {
            const volScalarField& alpha1 = iter1()();

            phaseModelTable::const_iterator iter2 = iter1;
            ++iter2;

            for (; iter2 != phaseModels_.end(); ++iter2)
            {
                const volScalarField& alpha2 = iter2()();

                stf +=
                    fvc::interpolate
                    (
                        surfaceTensionCoeff
                        (
                            phasePairKey(iter1()->name(), iter2()->name())
                        )
                    )
                * fvc::interpolate(K(alpha1, alpha2))*
                    (
                        fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                      - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                    );
            }
        }
    }

    return tstf;
}


Foam::tmp<Foam::volVectorField> Foam::phaseSystem::U() const
{
    auto tstf = tmp<volVectorField>::New
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector("U", dimVelocity, Zero)
    );

    auto& stf = tstf.ref();

    forAllConstIter(phaseModelTable, phaseModels_, iter1)
    {
        stf += iter1()()*iter1()->U();
    }

    return tstf;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::surfaceTensionCoeff(const phasePairKey& key) const
{
    return surfaceTensionModels_[key]->sigma();
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::coeffs
(
    const word& key
) const
{
    return 1.0/(phaseModels_[key]->thermo().rho());
}


void Foam::phaseSystem::addInterfacePorosity(fvVectorMatrix& UEqn)
{
    const scalarField& Vc = mesh_.V();
    scalarField& Udiag = UEqn.diag();

    forAllIter(phaseModelTable,  phaseModels_, iteri)
    {
        const phaseModel& phasei = iteri();

        phaseModelTable::iterator iterk = iteri;
        ++iterk;
        for
        (
            ;
            iterk != phaseModels_.end();
            ++iterk
        )
        {
            if (iteri()().name() != iterk()().name())
            {
                phaseModel& phasek = iterk()();

                // Phase i and k
                const phasePairKey keyik
                (
                    phasei.name(),
                    phasek.name(),
                    false
                );

                if (interfacePorousModelTable_.found(keyik))
                {
                    autoPtr<porousModel>& interfacePtr =
                        interfacePorousModelTable_[keyik];

                    Udiag += Vc*interfacePtr->S();
                }
            }
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::K
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    // Simple expression for curvature
    return -fvc::div(tnHatfv.ref() & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::nearInterface
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    return
    (
        pos(alpha1 - 0.1)*pos(0.9 - alpha1)
       *pos(alpha2 - 0.1)*pos(0.9 - alpha2)
    );
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::nearInterface() const
{
    auto tnI = tmp<volScalarField>::New
    (
        IOobject
        (
            "nearInterface",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    auto& nI = tnI.ref();

    forAllConstIter(phaseModelTable, phaseModels_, iter1)
    {
        const volScalarField& alpha1 = iter1()();

        phaseModelTable::const_iterator iter2 = iter1;
        ++iter2;

        for (; iter2 != phaseModels_.end(); ++iter2)
        {
            const volScalarField& alpha2 = iter2()();

            nI +=
            (
                pos(alpha1 - 0.1)*pos(0.9 - alpha1)
               *pos(alpha2 - 0.1)*pos(0.9 - alpha2)
            );
        }
    }

    return tnI;
}


Foam::tmp<Foam::surfaceVectorField> Foam::phaseSystem::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN);
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


bool Foam::phaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        return readOK;
    }

    return false;
}


// ************************************************************************* //

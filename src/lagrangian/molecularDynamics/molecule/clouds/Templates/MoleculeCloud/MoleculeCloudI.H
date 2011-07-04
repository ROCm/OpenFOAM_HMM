/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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

#include "constants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::evaluatePair
(
    MoleculeType& molI,
    MoleculeType& molJ
)
{
    const pairPotentialList& pairPot = pot_.pairPotentials();

    const pairPotential& electrostatic = pairPot.electrostatic();

    label idI = molI.id();

    label idJ = molJ.id();

    const typename MoleculeType::constantProperties& constPropI
    (
        constProps(idI)
    );

    const typename MoleculeType::constantProperties& constPropJ
    (
        constProps(idJ)
    );

    forAll(constPropI.pairPotSites(), pI)
    {
        label sI = constPropI.pairPotSites()[pI];

        label idsI = constPropI.sites()[sI].siteId();

        forAll(constPropJ.pairPotSites(), pJ)
        {
            label sJ = constPropJ.pairPotSites()[pJ];

            label idsJ = constPropJ.sites()[sJ].siteId();

            vector rsIsJ =
                molI.sitePositions()[sI] - molJ.sitePositions()[sJ];

            scalar rsIsJMagSq = magSqr(rsIsJ);

            if (pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
            {
                scalar rsIsJMag = mag(rsIsJ);

                vector fsIsJ =
                (rsIsJ/rsIsJMag)
               *pairPot.force(idsI, idsJ, rsIsJMag);

                molI.siteForces()[sI] += fsIsJ;

                molJ.siteForces()[sJ] += -fsIsJ;

                scalar potentialEnergy
                (
                    pairPot.energy(idsI, idsJ, rsIsJMag)
                );

                molI.potentialEnergy() += 0.5*potentialEnergy;

                molJ.potentialEnergy() += 0.5*potentialEnergy;

                vector rIJ = molI.position() - molJ.position();

                tensor virialContribution =
                (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;

                molI.rf() += virialContribution;

                molJ.rf() += virialContribution;
            }
        }
    }


    forAll(constPropI.electrostaticSites(), pI)
    {
        label sI = constPropI.electrostaticSites()[pI];

        forAll(constPropJ.electrostaticSites(), pJ)
        {
            label sJ = constPropJ.electrostaticSites()[pJ];

            vector rsIsJ =
                molI.sitePositions()[sI] - molJ.sitePositions()[sJ];

            scalar rsIsJMagSq = magSqr(rsIsJ);

            if (rsIsJMagSq <= electrostatic.rCutSqr())
            {
                scalar rsIsJMag = mag(rsIsJ);

                scalar chargeI = constPropI.sites()[sI].siteCharge();

                scalar chargeJ = constPropJ.sites()[sJ].siteCharge();

                vector fsIsJ =
                    (rsIsJ/rsIsJMag)
                   *chargeI*chargeJ*electrostatic.force(rsIsJMag);

                molI.siteForces()[sI] += fsIsJ;

                molJ.siteForces()[sJ] += -fsIsJ;

                scalar potentialEnergy =
                    chargeI*chargeJ
                   *electrostatic.energy(rsIsJMag);

                molI.potentialEnergy() += 0.5*potentialEnergy;

                molJ.potentialEnergy() += 0.5*potentialEnergy;

                vector rIJ = molI.position() - molJ.position();

                tensor virialContribution =
                    (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;

                molI.rf() += virialContribution;

                molJ.rf() += virialContribution;
            }
        }
    }
}


template<class MoleculeType>
bool Foam::MoleculeCloud<MoleculeType>::evaluatePotentialLimit
(
    MoleculeType& molI,
    MoleculeType& molJ
) const
{
    const pairPotentialList& pairPot = pot_.pairPotentials();

    const pairPotential& electrostatic = pairPot.electrostatic();

    label idI = molI.id();

    label idJ = molJ.id();

    const typename MoleculeType::constantProperties& constPropI
    (
        constProps(idI)
    );

    const typename MoleculeType::constantProperties& constPropJ
    (
        constProps(idJ)
    );

    forAll(constPropI.pairPotSites(), pI)
    {
        label sI = constPropI.pairPotSites()[pI];

        label idsI = constPropI.sites()[sI].siteId();

        forAll(constPropJ.pairPotSites(), pJ)
        {
            label sJ = constPropJ.pairPotSites()[pJ];

            label idsJ = constPropJ.sites()[sJ].siteId();

            vector rsIsJ =
                molI.sitePositions()[sI] - molJ.sitePositions()[sJ];

            scalar rsIsJMagSq = magSqr(rsIsJ);

            if (pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
            {
                scalar rsIsJMag = mag(rsIsJ);

                // Guard against pairPot.energy being evaluated
                // if rIJMag < SMALL. A floating point exception will
                // happen otherwise.

                if (rsIsJMag < SMALL)
                {
                    WarningIn
                    (
                        "MoleculeCloud<MoleculeType>::"
                        "removeHighEnergyOverlaps()"
                    )
                        << "Molecule site pair closer than "
                        << SMALL
                        << ": mag separation = " << rsIsJMag
                        << ". These may have been placed on top of each"
                        << " other by a rounding error in mdInitialise in"
                        << " parallel or a block filled with moleculess"
                        << " twice. Removing one of the molecules."
                        << endl;

                    return true;
                }

                // Guard against pairPot.energy being evaluated if rIJMag <
                // rMin.  A tabulation lookup error will occur otherwise.

                if (rsIsJMag < pairPot.rMin(idsI, idsJ))
                {
                    return true;
                }

                if
                (
                    mag(pairPot.energy(idsI, idsJ, rsIsJMag))
                  > pot_.potentialEnergyLimit()
                )
                {
                    return true;
                };
            }
        }
    }

    forAll(constPropI.electrostaticSites(), pI)
    {
        label sI = constPropI.electrostaticSites()[pI];

        forAll(constPropJ.electrostaticSites(), pJ)
        {
            label sJ = constPropJ.electrostaticSites()[pJ];

            vector rsIsJ =
                molI.sitePositions()[sI] - molJ.sitePositions()[sJ];

            scalar rsIsJMagSq = magSqr(rsIsJ);

            if (pairPot.rCutMaxSqr(rsIsJMagSq))
            {
                scalar rsIsJMag = mag(rsIsJ);

                // Guard against pairPot.energy being evaluated
                // if rIJMag < SMALL. A floating point exception will
                // happen otherwise.

                if (rsIsJMag < SMALL)
                {
                    WarningIn
                    (
                        "MoleculeCloud<MoleculeType>::"
                        "removeHighEnergyOverlaps()"
                    )
                        << "Molecule site pair closer than "
                        << SMALL
                        << ": mag separation = " << rsIsJMag
                        << ". These may have been placed on top of each"
                        << " other by a rounding error in mdInitialise in"
                        << " parallel or a block filled with molecules"
                        << " twice. Removing one of the molecules."
                        << endl;

                    return true;
                }

                if (rsIsJMag < electrostatic.rMin())
                {
                    return true;
                }

                scalar chargeI = constPropI.sites()[sI].siteCharge();

                scalar chargeJ = constPropJ.sites()[sJ].siteCharge();

                if
                (
                    mag(chargeI*chargeJ*electrostatic.energy(rsIsJMag))
                  > pot_.potentialEnergyLimit()
                )
                {
                    return true;
                };
            }
        }
    }

    return false;
}


template<class MoleculeType>
Foam::vector
Foam::MoleculeCloud<MoleculeType>::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return sqrt(physicoChemical::k.value()*temperature/mass)*vector
    (
        rndGen_.GaussNormal(),
        rndGen_.GaussNormal(),
        rndGen_.GaussNormal()
    );
}


template<class MoleculeType>
Foam::vector
Foam::MoleculeCloud<MoleculeType>::equipartitionAngularMomentum
(
    scalar temperature,
    const typename MoleculeType::constantProperties& cP
)
{
    scalar sqrtKbT = sqrt(physicoChemical::k.value()*temperature);

    if (cP.linearMolecule())
    {
        return sqrtKbT*vector
        (
            0.0,
            sqrt(cP.momentOfInertia().yy())*rndGen_.GaussNormal(),
            sqrt(cP.momentOfInertia().zz())*rndGen_.GaussNormal()
        );
    }
    else
    {
        return sqrtKbT*vector
        (
            sqrt(cP.momentOfInertia().xx())*rndGen_.GaussNormal(),
            sqrt(cP.momentOfInertia().yy())*rndGen_.GaussNormal(),
            sqrt(cP.momentOfInertia().zz())*rndGen_.GaussNormal()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MoleculeType>
const Foam::polyMesh& Foam::MoleculeCloud<MoleculeType>::mesh() const
{
    return mesh_;
}


template<class MoleculeType>
const Foam::potential& Foam::MoleculeCloud<MoleculeType>::pot() const
{
    return pot_;
}


template<class MoleculeType>
const Foam::List<Foam::DynamicList<MoleculeType*> >&
Foam::MoleculeCloud<MoleculeType>::cellOccupancy() const
{
    return cellOccupancy_;
}


template<class MoleculeType>
const Foam::InteractionLists<MoleculeType>&
Foam::MoleculeCloud<MoleculeType>::il() const
{
    return il_;
}


template<class MoleculeType>
const Foam::List<typename MoleculeType::constantProperties>
Foam::MoleculeCloud<MoleculeType>::constProps() const
{
    return constPropList_;
}


template<class MoleculeType>
const typename MoleculeType::constantProperties&
Foam::MoleculeCloud<MoleculeType>::constProps(label id) const
{
    return constPropList_[id];
}


template<class MoleculeType>
Foam::Random& Foam::MoleculeCloud<MoleculeType>::rndGen()
{
    return rndGen_;
}


// ************************************************************************* //

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

inline void Foam::polyatomicCloud::evaluatePair
(
    polyatomic& molI,
    polyatomic& molJ
)
{
    const pairPotentialList& pairPot = pot_.pairPotentials();

    const pairPotential& electrostatic = pairPot.electrostatic();

    label idI = molI.id();

    label idJ = molJ.id();

    const polyatomic::constantProperties& constPropI(constProps(idI));

    const polyatomic::constantProperties& constPropJ(constProps(idJ));

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


inline bool Foam::polyatomicCloud::evaluatePotentialLimit
(
    polyatomic& molI,
    polyatomic& molJ
) const
{
    const pairPotentialList& pairPot = pot_.pairPotentials();

    const pairPotential& electrostatic = pairPot.electrostatic();

    label idI = molI.id();

    label idJ = molJ.id();

    const polyatomic::constantProperties& constPropI(constProps(idI));

    const polyatomic::constantProperties& constPropJ(constProps(idJ));

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
                    WarningIn("polyatomicCloud::removeHighEnergyOverlaps()")
                        << "Molecule site pair closer than "
                        << SMALL
                        << ": mag separation = " << rsIsJMag
                        << ". These may have been placed on top of each"
                        << " other by a rounding error in mdInitialise in"
                        << " parallel or a block filled with polyatomics"
                        << " twice. Removing one of the polyatomics."
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
                    WarningIn("polyatomicCloud::removeHighEnergyOverlaps()")
                        << "Molecule site pair closer than "
                        << SMALL
                        << ": mag separation = " << rsIsJMag
                        << ". These may have been placed on top of each"
                        << " other by a rounding error in mdInitialise in"
                        << " parallel or a block filled with polyatomics"
                        << " twice. Removing one of the polyatomics."
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


inline Foam::vector Foam::polyatomicCloud::equipartitionLinearVelocity
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


inline Foam::vector Foam::polyatomicCloud::equipartitionAngularMomentum
(
    scalar temperature,
    const polyatomic::constantProperties& cP
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

inline const Foam::polyMesh& Foam::polyatomicCloud::mesh() const
{
    return mesh_;
}


inline const Foam::potential& Foam::polyatomicCloud::pot() const
{
    return pot_;
}


inline const Foam::List<Foam::DynamicList<Foam::polyatomic*> >&
    Foam::polyatomicCloud::cellOccupancy() const
{
    return cellOccupancy_;
}


inline const Foam::InteractionLists<Foam::polyatomic>&
    Foam::polyatomicCloud::il() const
{
    return il_;
}


inline const Foam::List<Foam::polyatomic::constantProperties>
    Foam::polyatomicCloud::constProps() const
{
    return constPropList_;
}


inline const Foam::polyatomic::constantProperties&
    Foam::polyatomicCloud::constProps(label id) const
{
    return constPropList_[id];
}


inline Foam::Random& Foam::polyatomicCloud::rndGen()
{
    return rndGen_;
}


// ************************************************************************* //

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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::moleculeCloud::evaluatePair
(
    molecule* molI,
    molecule* molJ
)
{
    const pairPotentialList& pairPot(pot_.pairPotentials());

    const electrostaticPotential& electrostatic(pot_.electrostatic());

    label idI = molI->id();

    label idJ = molJ->id();

    List<label> siteIdsI = constProps(idI).siteIds();

    List<label> siteIdsJ = constProps(idJ).siteIds();

    List<bool> pairPotentialSitesI = constProps(idI).pairPotentialSites();

    List<bool> electrostaticSitesI = constProps(idI).electrostaticSites();

    List<bool> pairPotentialSitesJ = constProps(idJ).pairPotentialSites();

    List<bool> electrostaticSitesJ = constProps(idJ).electrostaticSites();

    forAll(siteIdsI, sI)
    {
        label idsI(siteIdsI[sI]);

        forAll(siteIdsJ, sJ)
        {
            label idsJ(siteIdsJ[sJ]);

            if (pairPotentialSitesI[sI] && pairPotentialSitesJ[sJ])
            {
                Info<< "Pair potential " << idsI << ' ' << idsJ << endl;

                vector rsIsJ =
                molI->sitePositions()[sI] - molJ->sitePositions()[sJ];

                scalar rsIsJMagSq = magSqr(rsIsJ);

                if(pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                {
                    scalar rsIsJMag = mag(rsIsJ);

                    vector fsIsJ = (rsIsJ/rsIsJMag)
                    *pairPot.force(idsI, idsJ, rsIsJMag);

                    molI->siteForces()[sI] += fsIsJ;

                    molJ->siteForces()[sJ] += -fsIsJ;

                    scalar potentialEnergy
                    (
                        pairPot.energy(idsI, idsJ, rsIsJMag)
                    );

                    molI->potentialEnergy() += 0.5*potentialEnergy;

                    molJ->potentialEnergy() += 0.5*potentialEnergy;

                    // What to do here?

                    // molI->rf() += rIJ * fIJ;

                    // molJ->rf() += rIJ * fIJ;
                }
            }

            if (electrostaticSitesI[sI] && electrostaticSitesJ[sJ])
            {
                Info<< "Electrostatic " << idsI << ' ' << idsJ << endl;

                vector rsIsJ =
                molI->sitePositions()[sI] - molJ->sitePositions()[sJ];

                scalar rsIsJMag = mag(rsIsJ);

                vector fsIsJ = (rsIsJ/rsIsJMag)
                *constProps(idI).siteCharges()[sI]
                *constProps(idJ).siteCharges()[sJ]
                *electrostatic.force(rsIsJMag);

                molI->siteForces()[sI] += fsIsJ;

                molJ->siteForces()[sJ] += -fsIsJ;

                scalar potentialEnergy(electrostatic.energy(rsIsJMag));

                molI->potentialEnergy() += 0.5*potentialEnergy;

                molJ->potentialEnergy() += 0.5*potentialEnergy;

                // What to do here?

                // molI->rf() += rIJ * fIJ;

                // molJ->rf() += rIJ * fIJ;
            }
        }
    }
}


// inline void Foam::moleculeCloud::evaluatePair
// (
//     molecule* molI,
//     molecule* molJ
// )
// {
//     const pairPotentialList& pairPot(pot_.pairPotentials());

//     label idI = molI->id();

//     label idJ = molJ->id();

//     vector rIJ = molI->position() - molJ->position();

//     scalar rIJMagSq = magSqr(rIJ);

//     if(pairPot.rCutSqr(idI, idJ, rIJMagSq))
//     {
//         scalar rIJMag = mag(rIJ);

//         vector fIJ = (rIJ/rIJMag)*pairPot.force(idI, idJ, rIJMag);

//         scalar massI = constProps(idI).mass();

//         scalar massJ = constProps(idJ).mass();

//         // Acceleration increment for molI
//         molI->a() += fIJ/massI;

//         // Acceleration increment for molJ
//         molJ->a() += -fIJ/massJ;

//         scalar potentialEnergy
//         (
//             pairPot.energy(idI, idJ, rIJMag)
//         );

//         molI->potentialEnergy() += 0.5*potentialEnergy;

//         molJ->potentialEnergy() += 0.5*potentialEnergy;

//         molI->rf() += rIJ * fIJ;

//         molJ->rf() += rIJ * fIJ;
//     }
// }

inline void Foam::moleculeCloud::evaluatePair
(
    molecule* molReal,
    referredMolecule* molRef
)
{
    const pairPotentialList& pairPot(pot_.pairPotentials());

    label idReal = molReal->id();

    label idRef = molRef->id();

    vector rRealRef = molReal->position() - molRef->position();

    scalar rRealRefMagSq = magSqr(rRealRef);

    if (pairPot.rCutSqr(idReal, idRef, rRealRefMagSq))
    {
        scalar rRealRefMag = mag(rRealRef);

        vector fRealRef = (rRealRef/rRealRefMag)
            *pairPot.force(idReal, idRef, rRealRefMag);

        scalar massReal = constProps(idReal).mass();

        // Acceleration increment for molReal
        molReal->a() += fRealRef/massReal;

        // Adding a contribution of 1/2 of the potential energy
        // from this interaction
        molReal->potentialEnergy() +=
            0.5*pairPot.energy(idReal, idRef, rRealRefMag);

        molReal->rf() += rRealRef * fRealRef;
    }
}


inline bool Foam::moleculeCloud::evaluatePotentialLimit
(
    label idI,
    label idJ,
    scalar rIJMagSq
) const
{
    Info << "Move Foam::moleculeCloud::evaluatePotentialLimit "
        << "(label idI, label idJ, scalar rIJMagSq) const "
        << "to pairPotential class." << endl;

    const pairPotentialList& pairPot(pot_.pairPotentials());

    bool remove = false;

    if (pairPot.rCutSqr(idI, idJ, rIJMagSq))
    {
        scalar rIJMag = sqrt(rIJMagSq);

        // Guard against pairPot.energy being evaluated
        // if rIJMag < SMALL. A floating point exception will
        // happen otherwise.

        if (rIJMag < SMALL)
        {
            WarningIn("moleculeCloud::removeHighEnergyOverlaps()")
            << "Molecule pair are closer than "
                << SMALL
                << ": mag separation = " << rIJMag
                << ". These may have been placed on top of each"
                << " other by a rounding error in molConfig in parallel"
                << " or a block filled with molecules twice."
                << " Removing one of the molecules."
                << endl;

            remove = true;
        }

        // Guard against pairPot.energy being evaluated if rIJMag < rMin.
        // A tabulation lookup error will occur otherwise.

        if (rIJMag < pairPot.rMin(idI, idJ))
        {
            remove = true;
        }

        if (!remove)
        {
            if
            (
                pairPot.energy(idI, idJ, rIJMag) > pot_.potentialEnergyLimit()
            )
            {
                remove = true;
            }
        }
    }

    return remove;
}


inline bool Foam::moleculeCloud::evaluatePotentialLimit
(
    molecule* molI,
    molecule* molJ
) const
{
    label idI = molI->id();

    label idJ = molJ->id();

    vector rIJ = molI->position() - molJ->position();

    scalar rIJMagSq = magSqr(rIJ);

    return evaluatePotentialLimit(idI, idJ, rIJMagSq);
}


inline bool Foam::moleculeCloud::evaluatePotentialLimit
(
    molecule* molReal,
    referredMolecule* molRef
) const
{

    label idReal = molReal->id();

    label idRef = molRef->id();

    vector rRealRef = molReal->position() - molRef->position();

    scalar rRealRefMagSq = magSqr(rRealRef);

    return evaluatePotentialLimit(idReal, idRef, rRealRefMagSq);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::polyMesh& Foam::moleculeCloud::mesh() const
{
    return mesh_;
}


inline const Foam::potential& Foam::moleculeCloud::pot() const
{
    return pot_;
}


inline const Foam::List<Foam::DynamicList<Foam::molecule*> >&
    Foam::moleculeCloud::cellOccupancy() const
{
    return cellOccupancy_;
}


inline const Foam::interactionLists&
    Foam::moleculeCloud::il() const
{
    return il_;
}


inline const Foam::List<Foam::molecule::constantProperties>
    Foam::moleculeCloud::constProps() const
{
    return constPropList_;
}


inline const Foam::molecule::constantProperties&
    Foam::moleculeCloud::constProps(label id) const
{
    return constPropList_[id];
}


// ************************************************************************* //

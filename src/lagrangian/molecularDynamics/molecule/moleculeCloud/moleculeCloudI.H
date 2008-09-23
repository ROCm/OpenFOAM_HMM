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

    idI = molI->id();

    idJ = molJ->id();

    rIJ = molI->position() - molJ->position();

    rIJMagSq = magSqr(rIJ);

    if(pairPot.rCutSqr(idI, idJ, rIJMagSq))
    {
        rIJMag = mag(rIJ);

        fIJ = (rIJ/rIJMag)*pairPot.force(idI, idJ, rIJMag);

        // Acceleration increment for molI
        molI->A() += fIJ/(molI->mass());

        // Acceleration increment for molJ
        molJ->A() += -fIJ/(molJ->mass());

        scalar potentialEnergy
        (
            pairPot.energy(idI, idJ, rIJMag)
        );

        molI->potentialEnergy() += 0.5*potentialEnergy;

        molJ->potentialEnergy() += 0.5*potentialEnergy;

        molI->rf() += rIJ * fIJ;

        molJ->rf() += rIJ * fIJ;
    }
}

inline void Foam::moleculeCloud::evaluatePair
(
    molecule* molReal,
    referredMolecule* molRef
)
{
    const pairPotentialList& pairPot(pot_.pairPotentials());

    idReal = molReal->id();

    idRef = molRef->id();

    rRealRef = molReal->position() - molRef->position();

    rRealRefMagSq = magSqr(rRealRef);

    if (pairPot.rCutSqr(idReal, idRef, rRealRefMagSq))
    {
        rRealRefMag = mag(rRealRef);

        fRealRef = (rRealRef/rRealRefMag)
            *pairPot.force(idReal, idRef, rRealRefMag);

        // Acceleration increment for molReal
        molReal->A() += fRealRef/(molReal->mass());

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
        << "to pairPotential class."

    const pairPotentialList& pairPot(pot_.pairPotentials());

    bool remove = false;

    if (pairPot.rCutSqr(idI, idJ, rIJMagSq))
    {
        rIJMag = mag(rIJ);

        // Guard against pairPot.energy being evaluated
        // if rIJMag < SMALL. A floating point exception will
        // happen otherwise.

        if (rIJMag < SMALL)
        {
            WarningIn("moleculeCloud::removeHighEnergyOverlaps()")
            << "Real molecule pair "
                << " idI = " << idI
                << " at position " << molI->position()
                << " idJ = " << idJ
                << " at position " << molJ->position()
                << " are closer than " << SMALL
                << ": mag separation = " << rIJMag
                << ". These may have been placed on top of each"
                << " other by a rounding error in molConfig in parallel"
                << " or a block filled with molecules twice."
                << " Removing one of the molecules."
                << endl;

            remove = true;
        }

        // Guard against pairPot.energy being
        // evaluated if rIJMag < rMin. A tabulation lookup
        // error will occur otherwise.

        if (rIJMag < pairPot.rMin(idI, idJ))
        {
            remove = true;
        }

        if (!remove)
        {
            if
            (
                pairPot.energy(idI, idJ, rIJMag)
              > potentialEnergyLimit_
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
    idI = molI->id();

    idJ = molJ->id();

    rIJ = molI->position() - molJ->position();

    rIJMagSq = magSqr(rIJ);

    return evaluatePotentialLimit(idI, idJ, rIJMagSq);
}


inline bool Foam::moleculeCloud::evaluatePotentialLimit
(
    molecule* molReal,
    referredMolecule* molRef
) const
{

    idReal = molReal->id();

    idRef = molRef->id();

    rRealRef = molReal->position() - molRef->position();

    rRealRefMagSq = magSqr(rRealRef);

    return evaluatePotentialRefimit(idReal, idRef, rRealRefMagSq);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::polyMesh& moleculeCloud::mesh() const
{
    return mesh_;
}


inline const Foam::potential& moleculeCloud::potential() const
{
    return pot_;
}


inline const Foam::List<DynamicList<molecule*> >&
    Foam::moleculeCloud::cellOccupancy() const
{
    return cellOccupancy_;
}


inline const Foam::interactionLists&
    Foam::moleculeCloud::interactionLists() const
{
    return il_;
}


inline const Foam::List<molecule::constantProperties> constProps() const
{
    return constPropList_;
}


inline const molecule::constantProperties& constProps(label id) const
{
    return constPropList_[id];
}


// ************************************************************************* //

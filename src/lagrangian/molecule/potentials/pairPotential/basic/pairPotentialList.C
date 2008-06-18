/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "pairPotentialList.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pairPotentialList::pairPotentialList()
:
    List<pairPotential> ()
{}


Foam::pairPotentialList::pairPotentialList
(
    const label nIds
)
:
    List<pairPotential> ((nIds * (nIds + 1))/2),
    nIds_(nIds)
{}


Foam::pairPotentialList::pairPotentialList
(
    const List<pairPotential>& pairPotentials,
    const label nIds
)
:
    List<pairPotential> (pairPotentials),
    nIds_(nIds)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pairPotentialList::~pairPotentialList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pairPotentialList::setSizeAndNIds (const label nIds)
{
    nIds_ = nIds;
    setSize((nIds * (nIds + 1))/2);
}


void pairPotentialList::addPotential
(
    const label a,
    const label b,
    const pairPotential& pot
)
{
    if (pairPotentialIndex (a,b) > size()-1)
    {
        FatalErrorIn
        (
            "Foam::pairPotentialList::addPotential "
            "(const label a, const label b, const pairPotential& pot)"
        )<< "Attempting to add a pairPotential with too high an index"
         << nl
         << "Check if the pairPotentialList has been constructed "
         << "with the number of ids to expect"
         << nl << abort(FatalError);
    }
    else
    {
        (*this)[pairPotentialIndex (a,b)] = pot;
    }
}


const pairPotential& pairPotentialList::pairPotentialFunction
(
    const label a,
    const label b
) const
{
    return (*this)[pairPotentialIndex (a,b)];
}


bool pairPotentialList::rCutSqr
(
    const label a,
    const label b,
    const scalar rIJMagSqr
) const
{
    if(rIJMagSqr <= rCutSqr (a,b))
    {
        return true;
    }
    else
    {
        return false;
    }
}


scalar pairPotentialList::rCutSqr
(
    const label a,
    const label b
) const
{
    return (*this)[pairPotentialIndex (a,b)].rCutSqr();
}


scalar pairPotentialList::rCut
(
    const label a,
    const label b
) const
{
    return (*this)[pairPotentialIndex (a,b)].rCut();
}


scalar pairPotentialList::force
(
    const label a,
    const label b,
    const scalar rIJMag
) const
{
    scalar f = (*this)[pairPotentialIndex (a,b)].forceLookup(rIJMag);

    return f;
}


scalar pairPotentialList::energy
(
    const label a,
    const label b,
    const scalar rIJMag
) const
{
    scalar e = (*this)[pairPotentialIndex (a,b)].energyLookup(rIJMag);

    return e;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

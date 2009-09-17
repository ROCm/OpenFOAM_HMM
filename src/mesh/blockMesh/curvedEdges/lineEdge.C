/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "error.H"
#include "lineEdge.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lineEdge, 0);
    addToRunTimeSelectionTable(curvedEdge, lineEdge, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lineEdge::lineEdge
(
    const pointField& points,
    const label start,
    const label end
)
:
    curvedEdge(points, start, end),
    startPoint_(points_[start_]),
    direction_(points_[end_] - points_[start_])
{}


Foam::lineEdge::lineEdge(const pointField& points, Istream& is)
:
    curvedEdge(points, is),
    startPoint_(points_[start_]),
    direction_(points_[end_] - points_[start_])
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::lineEdge::position(const scalar lambda) const
{
    if (lambda < 0 || lambda > 1)
    {
        FatalErrorIn("lineEdge::position(const scalar)")
            << "Parameter out of range, lambda = " << lambda
            << abort(FatalError);
    }

    return startPoint_ + lambda*direction_;
}


Foam::scalar Foam::lineEdge::length() const
{
    return mag(direction_);
}


// ************************************************************************* //

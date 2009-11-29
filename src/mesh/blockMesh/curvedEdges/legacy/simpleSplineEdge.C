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

#include "simpleSplineEdge.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleSplineEdge, 0);
    addToRunTimeSelectionTable(curvedEdge, simpleSplineEdge, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleSplineEdge::simpleSplineEdge
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& otherknots
)
:
    curvedEdge(points, start, end),
    BSpline(appendEndPoints(points, start, end, otherknots))
{}


Foam::simpleSplineEdge::simpleSplineEdge
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& otherknots,
    const vector& fstend,
    const vector& sndend
)
:
    curvedEdge(points, start, end),
    BSpline(appendEndPoints(points, start, end, otherknots), fstend, sndend)
{}


Foam::simpleSplineEdge::simpleSplineEdge(const pointField& points, Istream& is)
:
    curvedEdge(points, is),
    BSpline(appendEndPoints(points, start_, end_, pointField(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::simpleSplineEdge::position(const scalar mu) const
{
    return BSpline::position(mu);
}


Foam::scalar Foam::simpleSplineEdge::length() const
{
    notImplemented("simpleSplineEdge::length() const");
    return 1.0;
}


// ************************************************************************* //

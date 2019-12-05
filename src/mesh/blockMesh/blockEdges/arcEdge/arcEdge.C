/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "arcEdge.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{
    defineTypeNameAndDebug(arcEdge, 0);
    addToRunTimeSelectionTable(blockEdge, arcEdge, Istream);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::coordSystem::cylindrical Foam::blockEdges::arcEdge::calcAngle()
{
    const vector a = p2_ - p1_;
    const vector b = p3_ - p1_;

    // Find centre of arcEdge
    const scalar asqr = a & a;
    const scalar bsqr = b & b;
    const scalar adotb = a & b;

    const scalar denom = asqr*bsqr - adotb*adotb;

    if (mag(denom) < VSMALL)
    {
        FatalErrorInFunction
            << denom
            << abort(FatalError);
    }

    const scalar fact = 0.5*(bsqr - adotb)/denom;

    point centre = 0.5*a + fact*((a ^ b) ^ a);

    centre += p1_;

    // Find position vectors w.r.t. the arcEdge centre
    const vector r1(p1_ - centre);
    const vector r2(p2_ - centre);
    const vector r3(p3_ - centre);

    // Find angle (in degrees)
    angle_ = radToDeg(acos((r3 & r1)/(mag(r3) * mag(r1))));

    // Check if the vectors define an exterior or an interior arcEdge
    if (((r1 ^ r2) & (r1 ^ r3)) < 0.0)
    {
        angle_ = 360.0 - angle_;
    }

    vector arcAxis;

    if (angle_ <= 180.0)
    {
        arcAxis = r1 ^ r3;

        if (mag(arcAxis)/(mag(r1)*mag(r3)) < 0.001)
        {
            arcAxis = r1 ^ r2;
        }
    }
    else
    {
        arcAxis = r3 ^ r1;
    }

    radius_ = mag(r3);

    // The corresponding local cylindrical coordinate system (radians)
    return coordSystem::cylindrical("arc", centre, arcAxis, r1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockEdges::arcEdge::arcEdge
(
    const pointField& points,
    const label start,
    const label end,
    const point& pMid
)
:
    blockEdge(points, start, end),
    p1_(points_[start_]),
    p2_(pMid),
    p3_(points_[end_]),
    cs_(calcAngle())
{}


Foam::blockEdges::arcEdge::arcEdge
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces& geometry,
    const pointField& points,
    Istream& is
)
:
    blockEdge(dict, index, points, is),
    p1_(points_[start_]),
    p2_(is),
    p3_(points_[end_]),
    cs_(calcAngle())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::blockEdges::arcEdge::position(const scalar lambda) const
{
    if (lambda < -SMALL || lambda > 1 + SMALL)
    {
        FatalErrorInFunction
            << "Parameter out of range, lambda = " << lambda
            << abort(FatalError);
    }

    if (lambda < SMALL)
    {
        return p1_;
    }
    else if (lambda > 1 - SMALL)
    {
        return p3_;
    }

    // The angle is degrees, the coordinate system in radians
    return cs_.globalPosition(vector(radius_, degToRad(lambda*angle_), 0));
}


Foam::scalar Foam::blockEdges::arcEdge::length() const
{
    return degToRad(radius_*angle_);
}


// ************************************************************************* //

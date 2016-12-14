/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "bezier.H"
#include "SubList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{
    defineTypeNameAndDebug(bezier, 0);
    addToRunTimeSelectionTable(blockEdge, bezier, Istream);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockEdges::bezier::bezier
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& control
)
:
    blockEdge(points, start, end),
    control_(control)
{}


Foam::blockEdges::bezier::bezier
(
    const dictionary& dict,
    const label index,
    const searchableSurfaces& geometry,
    const pointField& points,
    Istream& is
)
:
    blockEdge(dict, index, points, is),
    control_(appendEndPoints(points, start_, end_, pointField(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockEdges::bezier::~bezier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::blockEdges::bezier::position(const scalar lambda) const
{
    pointField working(control_);

    label nWorking(working.size());

    forAll(working, workingI)
    {
        --nWorking;

        SubList<point>(working, nWorking) =
            (1 - lambda)*SubList<point>(working, nWorking)
          + lambda*SubList<point>(working, nWorking, 1);
    }

    return working[0];
}


Foam::scalar Foam::blockEdges::bezier::length() const
{
    NotImplemented;
    return 1.0;
}


// ************************************************************************* //

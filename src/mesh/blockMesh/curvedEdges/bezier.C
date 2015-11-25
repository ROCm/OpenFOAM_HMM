/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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
    defineTypeNameAndDebug(bezier, 0);
    addToRunTimeSelectionTable(curvedEdge, bezier, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bezier::bezier
(
    const pointField& ps,
    const label start,
    const label end,
    const pointField& control
)
:
    curvedEdge(ps, start, end),
    control_(control)
{}


Foam::bezier::bezier(const pointField& ps, Istream& is)
:
    curvedEdge(ps, is),
    control_(appendEndPoints(ps, start_, end_, pointField(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bezier::~bezier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::bezier::position(const scalar lambda) const
{
    pointField working(control_);

    label nWorking(working.size());

    forAll(working, workingI)
    {
        -- nWorking;

        SubList<point>(working, nWorking).assign
        (
            (1 - lambda)*SubList<point>(working, nWorking)
          + lambda*SubList<point>(working, nWorking, 1)
        );
    }

    return working[0];
}


Foam::scalar Foam::bezier::length() const
{
    notImplemented("bezier::length() const");
    return 1.0;
}


// ************************************************************************* //

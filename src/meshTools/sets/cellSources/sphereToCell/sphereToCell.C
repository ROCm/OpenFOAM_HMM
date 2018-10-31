/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "sphereToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sphereToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, sphereToCell, word);
    addToRunTimeSelectionTable(topoSetSource, sphereToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, sphereToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, sphereToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        sphereToCell,
        word,
        sphere
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        sphereToCell,
        istream,
        sphere
    );
}


Foam::topoSetSource::addToUsageTable Foam::sphereToCell::usage_
(
    sphereToCell::typeName,
    "\n    Usage: sphereToCell (centreX centreY centreZ) radius\n\n"
    "    Select all cells with cellCentre within bounding sphere\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sphereToCell::combine(topoSet& set, const bool add) const
{
    const pointField& ctrs = mesh_.cellCentres();

    const scalar rad2 = radius_*radius_;

    forAll(ctrs, celli)
    {
        if (magSqr(ctrs[celli] - origin_) <= rad2)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sphereToCell::sphereToCell
(
    const polyMesh& mesh,
    const point& origin,
    const scalar radius
)
:
    topoSetCellSource(mesh),
    origin_(origin),
    radius_(radius)
{}


Foam::sphereToCell::sphereToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sphereToCell
    (
        mesh,
        dict.getCompat<vector>("origin", {{"centre", 1806}}),
        dict.get<scalar>("radius")
    )
{}


Foam::sphereToCell::sphereToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    origin_(checkIs(is)),
    radius_(readScalar(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sphereToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells within a sphere with centre = "
                << origin_ << " and radius = " << radius_ << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing cells within a sphere with centre = "
                << origin_ << " and radius = " << radius_ << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //

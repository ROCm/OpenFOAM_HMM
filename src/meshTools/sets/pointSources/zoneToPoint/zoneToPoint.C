/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "zoneToPoint.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneToPoint, 0);
    addToRunTimeSelectionTable(topoSetSource, zoneToPoint, word);
    addToRunTimeSelectionTable(topoSetSource, zoneToPoint, istream);
    addToRunTimeSelectionTable(topoSetPointSource, zoneToPoint, word);
    addToRunTimeSelectionTable(topoSetPointSource, zoneToPoint, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetPointSource,
        zoneToPoint,
        word,
        zone
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetPointSource,
        zoneToPoint,
        istream,
        zone
    );
}


Foam::topoSetSource::addToUsageTable Foam::zoneToPoint::usage_
(
    zoneToPoint::typeName,
    "\n    Usage: zoneToPoint zone\n\n"
    "    Select all points in the pointZone."
    " Note:accepts wildcards for zone.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneToPoint::combine(topoSet& set, const bool add) const
{
    bool hasMatched = false;

    for (const pointZone& zone : mesh_.pointZones())
    {
        if (selectedZones_.match(zone.name()))
        {
            hasMatched = true;

            const labelList& pointLabels = zone;

            if (verbose_)
            {
                Info<< "    Found matching zone " << zone.name()
                    << " with " << pointLabels.size() << " points." << endl;
            }

            for (const label pointi : pointLabels)
            {
                // Only do active points
                if (pointi >= 0 && pointi < mesh_.nPoints())
                {
                    addOrDelete(set, pointi, add);
                }
            }
        }
    }

    if (!hasMatched)
    {
        WarningInFunction
            << "Cannot find any pointZone matching "
            << flatOutput(selectedZones_) << nl
            << "Valid names: " << flatOutput(mesh_.pointZones().names())
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneToPoint::zoneToPoint
(
    const polyMesh& mesh,
    const wordRe& zoneName
)
:
    topoSetPointSource(mesh),
    selectedZones_(one(), zoneName)
{}


Foam::zoneToPoint::zoneToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetPointSource(mesh),
    selectedZones_()
{
    // Look for 'zones' and 'zone', but accept 'name' as well
    if (!dict.readIfPresent("zones", selectedZones_))
    {
        selectedZones_.resize(1);
        selectedZones_.first() =
            dict.getCompat<wordRe>("zone", {{"name", 1806}});
    }
}


Foam::zoneToPoint::zoneToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetPointSource(mesh),
    selectedZones_(one(), wordRe(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all points of point zones "
                << flatOutput(selectedZones_) << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all points of point zones "
                << flatOutput(selectedZones_) << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //

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

#include "zoneToPoint.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneToPoint, 0);
    addToRunTimeSelectionTable(topoSetSource, zoneToPoint, word);
    addToRunTimeSelectionTable(topoSetSource, zoneToPoint, istream);
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
        if (zoneName_.match(zone.name()))
        {
            hasMatched = true;

            const labelList& pointLabels = zone;

            Info<< "    Found matching zone " << zone.name()
                << " with " << pointLabels.size() << " points." << endl;

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
            << "Cannot find any pointZone named " << zoneName_ << nl
            << "Valid names: " << flatOutput(mesh_.pointZones().names())
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneToPoint::zoneToPoint
(
    const polyMesh& mesh,
    const word& zoneName
)
:
    topoSetSource(mesh),
    zoneName_(zoneName)
{}


Foam::zoneToPoint::zoneToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    zoneName_(dict.get<wordRe>("name"))
{}


Foam::zoneToPoint::zoneToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    zoneName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all points of pointZone " << zoneName_ << " ..."
            << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all points of pointZone " << zoneName_ << " ..."
            << endl;

        combine(set, false);
    }
}


// ************************************************************************* //

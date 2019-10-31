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

#include "zoneToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, zoneToCell, word);
    addToRunTimeSelectionTable(topoSetSource, zoneToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, zoneToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, zoneToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        zoneToCell,
        word,
        zone
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        zoneToCell,
        istream,
        zone
    );
}


Foam::topoSetSource::addToUsageTable Foam::zoneToCell::usage_
(
    zoneToCell::typeName,
    "\n    Usage: zoneToCell zone\n\n"
    "    Select all cells in the cellZone."
    " Note:accepts wildcards for zone.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneToCell::combine(topoSet& set, const bool add) const
{
    bool hasMatched = false;

    for (const cellZone& zone : mesh_.cellZones())
    {
        if (selectedZones_.match(zone.name()))
        {
            hasMatched = true;

            const labelList& cellLabels = zone;

            if (verbose_)
            {
                Info<< "    Found matching zone " << zone.name()
                    << " with " << cellLabels.size() << " cells." << endl;
            }

            for (const label celli : cellLabels)
            {
                // Only do active cells
                if (celli >= 0 && celli < mesh_.nCells())
                {
                    addOrDelete(set, celli, add);
                }
            }
        }
    }

    if (!hasMatched)
    {
        WarningInFunction
            << "Cannot find any cellZone matching "
            << flatOutput(selectedZones_) << nl
            << "Valid names: " << flatOutput(mesh_.cellZones().names())
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneToCell::zoneToCell
(
    const polyMesh& mesh,
    const wordRe& zoneName
)
:
    topoSetCellSource(mesh),
    selectedZones_(one(), zoneName)
{}


Foam::zoneToCell::zoneToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh),
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


Foam::zoneToCell::zoneToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    selectedZones_(one(), wordRe(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all cells of cell zones "
                << flatOutput(selectedZones_) << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all cells of cell zones "
                << flatOutput(selectedZones_) << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //

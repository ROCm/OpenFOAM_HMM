/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "faceZoneToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZoneToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, faceZoneToCell, word);
    addToRunTimeSelectionTable(topoSetSource, faceZoneToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, faceZoneToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, faceZoneToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::faceZoneToCell::usage_
(
    faceZoneToCell::typeName,
    "\n    Usage: faceZoneToCell zone master|slave\n\n"
    "    Select master or slave side of the faceZone."
    " Note:accepts wildcards for zone.\n\n"
);


const Foam::Enum
<
    Foam::faceZoneToCell::faceAction
>
Foam::faceZoneToCell::faceActionNames_
({
    { faceAction::MASTER, "master" },
    { faceAction::SLAVE, "slave" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceZoneToCell::combine(topoSet& set, const bool add) const
{
    bool hasMatched = false;

    for (const faceZone& zone : mesh_.faceZones())
    {
        if (selectedZones_.match(zone.name()))
        {
            hasMatched = true;

            const labelList& cellLabels =
            (
                option_ == MASTER
              ? zone.masterCells()
              : zone.slaveCells()
            );

            if (verbose_)
            {
                Info<< "    Found matching zone " << zone.name()
                    << " with " << cellLabels.size() << " cells on "
                    << faceActionNames_[option_] << " side" << endl;
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
            << "Cannot find any faceZone matching "
            << flatOutput(selectedZones_) << nl
            << "Valid names: " << flatOutput(mesh_.faceZones().names())
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZoneToCell::faceZoneToCell
(
    const polyMesh& mesh,
    const wordRe& zoneName,
    const faceAction option
)
:
    topoSetCellSource(mesh),
    selectedZones_(one(), zoneName),
    option_(option)
{}


Foam::faceZoneToCell::faceZoneToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh),
    selectedZones_(),
    option_(faceActionNames_.get("option", dict))
{
    // Look for 'zones' and 'zone', but accept 'name' as well
    if (!dict.readIfPresent("zones", selectedZones_))
    {
        selectedZones_.resize(1);
        selectedZones_.first() =
            dict.getCompat<wordRe>("zone", {{"name", 1806}});
    }
}


Foam::faceZoneToCell::faceZoneToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    selectedZones_(one(), wordRe(checkIs(is))),
    option_(faceActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceZoneToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all " << faceActionNames_[option_]
                << " cells of face zones "
                << flatOutput(selectedZones_) << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all " << faceActionNames_[option_]
                << " cells of face zones "
                << flatOutput(selectedZones_) << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //

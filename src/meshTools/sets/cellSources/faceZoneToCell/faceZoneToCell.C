/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    //TBD: { faceAction::BOTH, "attached" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceZoneToCell::combine
(
    topoSet& set,
    const labelUList& zoneIDs,
    const bool add,
    const bool verbosity
) const
{
    const label nZones = mesh_.faceZones().size();

    if (zoneIDs.empty() || !nZones)
    {
        return;  // Nothing to do
    }

    for (const label zonei : zoneIDs)
    {
        if (zonei < 0 || zonei >= nZones)
        {
            continue;
        }

        const auto& zone = mesh_.faceZones()[zonei];

        const labelList& cellLabels =
        (
            option_ == MASTER
          ? zone.masterCells()
          : zone.slaveCells()
        );

        if (verbosity)
        {
            Info<< "    Using matching zone " << zone.name()
                << " with " << cellLabels.size() << " cells on "
                << faceActionNames_[option_] << " side" << endl;
        }

        // NOTE could also handle both sides directly if required

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


void Foam::faceZoneToCell::combine(topoSet& set, const bool add) const
{
    if (zoneMatcher_.empty())
    {
        return;  // Nothing to do
    }

    const labelList matched(mesh_.faceZones().indices(zoneMatcher_));

    if (matched.empty())
    {
        WarningInFunction
            << "Cannot find any faceZone matching "
            << flatOutput(zoneMatcher_) << nl
            << "Valid names are " << flatOutput(mesh_.faceZones().names())
            << endl;

        return;  // Nothing to do
    }

    combine(set, matched, add, verbose_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZoneToCell::faceZoneToCell
(
    const polyMesh& mesh,
    const wordRes& zoneSelector,
    const faceAction option
)
:
    topoSetCellSource(mesh),
    zoneMatcher_(zoneSelector),
    option_(option)
{}


Foam::faceZoneToCell::faceZoneToCell
(
    const polyMesh& mesh,
    const wordRe& zoneName,
    const faceAction option
)
:
    topoSetCellSource(mesh),
    zoneMatcher_(one{}, zoneName),
    option_(option)
{}


Foam::faceZoneToCell::faceZoneToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh),
    zoneMatcher_(),
    option_(faceActionNames_.get("option", dict))
{
    // Look for 'zones' and 'zone', but accept 'name' as well
    if (!dict.readIfPresent("zones", zoneMatcher_))
    {
        zoneMatcher_.resize(1);
        zoneMatcher_.first() =
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
    zoneMatcher_(one{}, wordRe(checkIs(is))),
    option_(faceActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordRes& Foam::faceZoneToCell::zones() const noexcept
{
    return zoneMatcher_;
}


void Foam::faceZoneToCell::zones(const wordRes& zonesSelector)
{
    zoneMatcher_ = zonesSelector;
}


void Foam::faceZoneToCell::zones(const wordRe& zoneName)
{
    zoneMatcher_.resize(1);
    zoneMatcher_.first() = zoneName;
}


void Foam::faceZoneToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_ && !zoneMatcher_.empty())
        {
            Info<< "    Adding all " << faceActionNames_[option_]
                << " cells of face zones "
                << flatOutput(zoneMatcher_) << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_ && !zoneMatcher_.empty())
        {
            Info<< "    Removing all " << faceActionNames_[option_]
                << " cells of face zones "
                << flatOutput(zoneMatcher_) << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //

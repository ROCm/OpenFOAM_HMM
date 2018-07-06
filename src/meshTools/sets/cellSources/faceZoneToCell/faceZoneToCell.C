/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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
{
    { faceAction::MASTER, "master" },
    { faceAction::SLAVE, "slave" },
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceZoneToCell::combine(topoSet& set, const bool add) const
{
    bool hasMatched = false;

    for (const faceZone& zone : mesh_.faceZones())
    {
        if (zoneName_.match(zone.name()))
        {
            hasMatched = true;

            const labelList& cellLabels =
            (
                option_ == MASTER
              ? zone.masterCells()
              : zone.slaveCells()
            );

            Info<< "    Found matching zone " << zone.name()
                << " with " << cellLabels.size() << " cells on selected side."
                << endl;

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
            << "Cannot find any faceZone named " << zoneName_ << nl
            << "Valid names: " << flatOutput(mesh_.faceZones().names())
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZoneToCell::faceZoneToCell
(
    const polyMesh& mesh,
    const word& zoneName,
    const faceAction option
)
:
    topoSetSource(mesh),
    zoneName_(zoneName),
    option_(option)
{}


Foam::faceZoneToCell::faceZoneToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    zoneName_(dict.get<wordRe>("name")),
    option_(faceActionNames_.lookup("option", dict))
{}


Foam::faceZoneToCell::faceZoneToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    zoneName_(checkIs(is)),
    option_(faceActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceZoneToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all " << faceActionNames_[option_]
            << " cells of faceZone " << zoneName_ << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all " << faceActionNames_[option_]
            << " cells of faceZone " << zoneName_ << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //

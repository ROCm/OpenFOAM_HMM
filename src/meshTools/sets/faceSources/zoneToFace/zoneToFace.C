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

#include "zoneToFace.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, zoneToFace, word);
    addToRunTimeSelectionTable(topoSetSource, zoneToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, zoneToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, zoneToFace, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        zoneToFace,
        word,
        zone
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        zoneToFace,
        istream,
        zone
    );
}


Foam::topoSetSource::addToUsageTable Foam::zoneToFace::usage_
(
    zoneToFace::typeName,
    "\n    Usage: zoneToFace zone\n\n"
    "    Select all faces in the faceZone."
    " Note:accepts wildcards for zone.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneToFace::combine(topoSet& set, const bool add) const
{
    bool hasMatched = false;

    for (const faceZone& zone : mesh_.faceZones())
    {
        if (selectedZones_.match(zone.name()))
        {
            hasMatched = true;

            const labelList& faceLabels = zone;

            if (verbose_)
            {
                Info<< "    Found matching zone " << zone.name()
                    << " with " << faceLabels.size() << " faces." << endl;
            }

            for (const label facei : faceLabels)
            {
                // Only do active faces
                if (facei >= 0 && facei < mesh_.nFaces())
                {
                    addOrDelete(set, facei, add);
                }
            }
        }
    }

    if (!hasMatched)
    {
        WarningInFunction
            << "Cannot find any faceZone matching "
            << flatOutput(selectedZones_) << nl
            << "Valid names are " << flatOutput(mesh_.faceZones().names())
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneToFace::zoneToFace
(
    const polyMesh& mesh,
    const wordRe& zoneName
)
:
    topoSetFaceSource(mesh),
    selectedZones_(one(), zoneName)
{}


Foam::zoneToFace::zoneToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh),
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


Foam::zoneToFace::zoneToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceSource(mesh),
    selectedZones_(one(), wordRe(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all faces of face zones "
                << flatOutput(selectedZones_) << " ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all faces of face zones "
                << flatOutput(selectedZones_) << " ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //

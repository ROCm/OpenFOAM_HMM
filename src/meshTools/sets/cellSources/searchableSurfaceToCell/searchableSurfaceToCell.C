/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "searchableSurfaceToCell.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSurfaceToCell, 0);
    addToRunTimeSelectionTable
    (
        topoSetSource,
        searchableSurfaceToCell,
        word
    );
    addToRunTimeSelectionTable
    (
        topoSetCellSource,
        searchableSurfaceToCell,
        word
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetSource,
        searchableSurfaceToCell,
        word,
        surface
    );
}


Foam::topoSetSource::addToUsageTable Foam::searchableSurfaceToCell::usage_
(
    searchableSurfaceToCell::typeName,
    "\n    Usage: searchableSurfaceToCell surface\n\n"
    "    Select cells with centre enclosed by the surface"
    "\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableSurfaceToCell::combine(topoSet& set, const bool add) const
{
    if (!surf_)
    {
        return;
    }
    const searchableSurface& s = *surf_;

    // Add cells within the enclosing volumes

    const label len = mesh_.nCells();

    List<volumeType> volTypes;

    s.getVolumeType(mesh_.cellCentres(), volTypes);

    for (label celli=0; celli < len; ++celli)
    {
        if (volTypes[celli] == volumeType::INSIDE)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceToCell::searchableSurfaceToCell
(
    const word& surfaceType,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh),
    surf_
    (
        searchableSurface::New
        (
            surfaceType,
            IOobject
            (
                dict.lookupOrDefault("name", mesh.objectRegistry::db().name()),
                mesh.time().constant(), // Instance
                "triSurface",           // Local
                mesh.time(),            // Registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    )
{
    // Check/warn for non-enclosed
    if (surf_ && !surf_->hasVolumeType())
    {
        WarningInFunction
            << nl << "The surface '" << surf_->name() << "' of type '"
            << surf_->type() << "' appears to be unclosed ... ignoring"
            << nl << endl;

        surf_.clear();
    }
}


Foam::searchableSurfaceToCell::searchableSurfaceToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    searchableSurfaceToCell
    (
        dict.get<word>("surface"),
        mesh,
        dict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::searchableSurfaceToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (!surf_ || !surf_->hasVolumeType())
    {
        return;
    }

    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells enclosed by searchableSurface"
                << "..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing cells enclosed by searchableSurface"
                << "..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //

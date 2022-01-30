/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "faceToPoint.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceToPoint, 0);
    addToRunTimeSelectionTable(topoSetSource, faceToPoint, word);
    addToRunTimeSelectionTable(topoSetSource, faceToPoint, istream);
    addToRunTimeSelectionTable(topoSetPointSource, faceToPoint, word);
    addToRunTimeSelectionTable(topoSetPointSource, faceToPoint, istream);
}

Foam::topoSetSource::addToUsageTable Foam::faceToPoint::usage_
(
    faceToPoint::typeName,
    "\n    Usage: faceToPoint <faceSet> all\n\n"
    "    Select all points of faces in the faceSet\n\n"
);

const Foam::Enum
<
    Foam::faceToPoint::faceAction
>
Foam::faceToPoint::faceActionNames_
({
    { faceAction::ALL, "all" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceToPoint::combine
(
    topoSet& set,
    const bool add,
    const word& setName
) const
{
    // Load the set
    faceSet loadedSet(mesh_, setName);
    const labelHashSet& faceLabels = loadedSet;

    // Add all points from faces in loadedSet
    for (const label facei : faceLabels)
    {
        const face& f = mesh_.faces()[facei];

        addOrDelete(set, f, add);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceToPoint::faceToPoint
(
    const polyMesh& mesh,
    const word& setName,
    const faceAction option
)
:
    topoSetPointSource(mesh),
    names_(one{}, setName),
    option_(option)
{}


Foam::faceToPoint::faceToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetPointSource(mesh),
    names_(),
    option_(faceActionNames_.get("option", dict))
{
    // Look for 'sets' or 'set'
    if (!dict.readIfPresent("sets", names_))
    {
        names_.resize(1);
        dict.readEntry("set", names_.first());
    }
}


Foam::faceToPoint::faceToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetPointSource(mesh),
    names_(one{}, word(checkIs(is))),
    option_(faceActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding points from face in faceSet "
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            combine(set, true, setName);
        }
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing points from face in faceSet "
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            combine(set, false, setName);
        }
    }
}


// ************************************************************************* //

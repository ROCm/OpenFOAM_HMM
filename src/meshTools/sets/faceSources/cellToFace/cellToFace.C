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

#include "cellToFace.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "Time.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, cellToFace, word);
    addToRunTimeSelectionTable(topoSetSource, cellToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, cellToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, cellToFace, istream);
}


Foam::topoSetSource::addToUsageTable Foam::cellToFace::usage_
(
    cellToFace::typeName,
    "\n    Usage: cellToFace <cellSet> all|both\n\n"
    "    Select -all : all faces of cells in the cellSet\n"
    "           -both: faces where both neighbours are in the cellSet\n\n"
);

const Foam::Enum
<
    Foam::cellToFace::cellAction
>
Foam::cellToFace::cellActionNames_
({
    { cellAction::ALL, "all" },
    { cellAction::BOTH, "both" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellToFace::combine
(
    topoSet& set,
    const bool add,
    const word& setName
) const
{
    // Load the set
    if (!exists(mesh_.time().path()/topoSet::localPath(mesh_, setName)))
    {
        SeriousError<< "Cannot load set "
            << setName << endl;
    }

    cellSet loadedSet(mesh_, setName);
    const labelHashSet& cellLabels = loadedSet;

    if (option_ == ALL)
    {
        // Add all faces from cell
        for (const label celli : cellLabels)
        {
            const labelList& cFaces = mesh_.cells()[celli];

            addOrDelete(set, cFaces, add);
        }
    }
    else if (option_ == BOTH)
    {
        // Add all faces whose both neighbours are in set.

        const label nInt = mesh_.nInternalFaces();
        const labelList& own = mesh_.faceOwner();
        const labelList& nei = mesh_.faceNeighbour();
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();


        // Check all internal faces
        for (label facei = 0; facei < nInt; ++facei)
        {
            if (cellLabels.found(own[facei]) && cellLabels.found(nei[facei]))
            {
                addOrDelete(set, facei, add);
            }
        }


        // Get coupled cell status
        boolList neiInSet(mesh_.nFaces()-nInt, false);

        for (const polyPatch& pp : patches)
        {
            if (pp.coupled())
            {
                label facei = pp.start();
                forAll(pp, i)
                {
                    neiInSet[facei-nInt] = cellLabels.found(own[facei]);
                    ++facei;
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiInSet);


        // Check all boundary faces
        for (const polyPatch& pp : patches)
        {
            if (pp.coupled())
            {
                label facei = pp.start();
                forAll(pp, i)
                {
                    if (cellLabels.found(own[facei]) && neiInSet[facei-nInt])
                    {
                        addOrDelete(set, facei, add);
                    }
                    ++facei;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellToFace::cellToFace
(
    const polyMesh& mesh,
    const word& setName,
    const cellAction option
)
:
    topoSetFaceSource(mesh),
    names_(one{}, setName),
    option_(option)
{}


Foam::cellToFace::cellToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh),
    names_(),
    option_(cellActionNames_.get("option", dict))
{
    // Look for 'sets' or 'set'
    if (!dict.readIfPresent("sets", names_))
    {
        names_.resize(1);
        dict.readEntry("set", names_.first());
    }
}


Foam::cellToFace::cellToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceSource(mesh),
    names_(one{}, word(checkIs(is))),
    option_(cellActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding faces according to cellSet "
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
            Info<< "    Removing faces according to cellSet "
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            combine(set, false, setName);
        }
    }
}


// ************************************************************************* //

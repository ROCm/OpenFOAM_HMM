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

#include "pointToCell.H"
#include "polyMesh.H"
#include "pointSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, pointToCell, word);
    addToRunTimeSelectionTable(topoSetSource, pointToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::pointToCell::usage_
(
    pointToCell::typeName,
    "\n    Usage: pointToCell <pointSet> any|edge\n\n"
    "    Select all cells with any point ('any') or any edge ('edge')"
    " in the pointSet\n\n"
);

const Foam::Enum
<
    Foam::pointToCell::pointAction
>
Foam::pointToCell::pointActionNames_
{
    { pointAction::ANY, "any" },
    { pointAction::EDGE, "edge" },
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointToCell::combine(topoSet& set, const bool add) const
{
    // Load the set
    pointSet loadedSet(mesh_, setName_);

    const labelHashSet& pointLabels = loadedSet;

    // Handle any selection
    if (option_ == ANY)
    {
        for (const label pointi : pointLabels)
        {
            const labelList& pCells = mesh_.pointCells()[pointi];

            for (const label celli : pCells)
            {
                addOrDelete(set, celli, add);
            }
        }
    }
    else if (option_ == EDGE)
    {
        const faceList& faces = mesh_.faces();

        forAll(faces, facei)
        {
            const face& f = faces[facei];

            forAll(f, fp)
            {
                if
                (
                    pointLabels.found(f[fp])
                 && pointLabels.found(f.nextLabel(fp))
                )
                {
                    addOrDelete(set, mesh_.faceOwner()[facei], add);
                    if (mesh_.isInternalFace(facei))
                    {
                        addOrDelete(set, mesh_.faceNeighbour()[facei], add);
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    const word& setName,
    const pointAction option
)
:
    topoSetSource(mesh),
    setName_(setName),
    option_(option)
{}


Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.get<word>("set")),
    option_(pointActionNames_.lookup("option", dict))
{}


Foam::pointToCell::pointToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    option_(pointActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells according to pointSet " << setName_
            << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells according to pointSet " << setName_
            << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //

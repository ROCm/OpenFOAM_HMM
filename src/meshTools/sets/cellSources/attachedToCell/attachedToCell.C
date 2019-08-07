/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

#include "attachedToCell.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(attachedToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, attachedToCell, word);
    addToRunTimeSelectionTable(topoSetSource, attachedToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, attachedToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, attachedToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        attachedToCell,
        word,
        attached
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        attachedToCell,
        istream,
        attached
    );
}


Foam::topoSetSource::addToUsageTable Foam::attachedToCell::usage_
(
    attachedToCell::typeName,
    "\n    Usage: attachedToCell\n\n"
    "    Select attached cells\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::attachedToCell::combine(topoSet& set, const bool add) const
{
    const cellList& cells = mesh_.cells();
    const labelList& faceOwn = mesh_.faceOwner();
    const labelList& faceNei = mesh_.faceNeighbour();
    const label nIntFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    boolList isCoupled(mesh_.nBoundaryFaces(), false);

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            label facei = pp.start();
            forAll(pp, i)
            {
                isCoupled[facei-nIntFaces] = true;
                ++facei;
            }
        }
    }

    // The starting set of cells
    bitSet current(cells.size());

    for (const label celli : set)
    {
        current.set(celli);
    }

    // The perimeter faces of the cell set
    bitSet outsideFaces(mesh().nFaces());

    // Get coupled cell status
    boolList neiInSet(mesh_.nFaces()-nIntFaces, false);

    bitSet updates(current);

    for (label stepi = 0; stepi < steps_; ++stepi)
    {
        // Mark up all perimeter faces

        // Don't attempt extra efficiency with caching old values,
        // this make the parallel transfer too troublesome

        outsideFaces.reset();
        for (const label celli : current)
        {
            for (const label facei : cells[celli])
            {
                outsideFaces.flip(facei);
            }
        }

        // Coupled faces
        neiInSet = false;

        for
        (
            bitSet::const_iterator iter = outsideFaces.cbegin(nIntFaces);
            iter != outsideFaces.cend();
            ++iter
        )
        {
            const label facei = *iter;

            if (isCoupled[facei-nIntFaces])
            {
                neiInSet[facei-nIntFaces] = true;
            }
        }

        syncTools::swapBoundaryFaceList(mesh_, neiInSet);

        forAll(neiInSet, bfacei)
        {
            if (neiInSet[bfacei])
            {
                outsideFaces.flip(bfacei+nIntFaces);
            }
        }

        updates.reset();

        for (const label facei : outsideFaces)
        {
            updates.set(faceOwn[facei]);
            if (mesh_.isInternalFace(facei))
            {
                updates.set(faceNei[facei]);
            }
        }

        if (add)
        {
            // Face information retained for next loop

            // Changed cells
            updates -= current;

            if (verbose_)
            {
                Info<< "    Grow " << current.count()
                    << " by " << updates.count() << endl;
            }

            current |= updates;
        }
        else
        {
            // Changed cells
            updates &= current;

            if (verbose_)
            {
                Info<< "    Shrink " << current.count()
                    << " by " << updates.count() << endl;
            }

            current -= updates;
        }

        for (const label celli: updates)
        {
            addOrDelete(set, celli, add);
        }

        if (updates.none())
        {
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::attachedToCell::attachedToCell
(
    const polyMesh& mesh,
    const label steps
)
:
    topoSetCellSource(mesh),
    steps_(max(steps, 1))
{}


Foam::attachedToCell::attachedToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    attachedToCell(mesh, dict.getOrDefault<label>("steps", 1))
{}


Foam::attachedToCell::attachedToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    attachedToCell(mesh, 1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::attachedToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Cannot create new of attached" << endl;
        }

        set.clear();
    }
    else if (action == topoSetSource::ADD)
    {
        if (verbose_)
        {
            Info<< "    Adding cells attached to current set, using "
                << steps_ << " step ..." << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing external cells of current set, using "
                << steps_ << " step ..." << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //

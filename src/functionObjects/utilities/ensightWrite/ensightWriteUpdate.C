/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "ensightWrite.H"
#include "dictionary.H"
#include "cellBitSet.H"
#include "topoSetCellSource.H"

// * * * * * * * * * * * * * * Local Data Members  * * * * * * * * * * * * * //

namespace Foam
{
    // A limited selection of actions
    const Enum<topoSetSource::setAction> actionNames
    ({
        { topoSetSource::NEW, "use" },  // Reuse NEW for "use" action name
        { topoSetSource::ADD, "add" },
        { topoSetSource::SUBTRACT, "subtract" },
        { topoSetSource::SUBSET, "subset" },
        { topoSetSource::INVERT, "invert" },
    });
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::ensightWrite::updateSubset
(
    fvMeshSubset& subsetter
) const
{
    if (selection_.empty())
    {
        return false;
    }

    const fvMesh& mesh = subsetter.baseMesh();

    // Start with all cells unselected
    cellBitSet cellsToSelect(mesh, false);

    // Execute all actions
    for (const entry& dEntry : selection_)
    {
        if (!dEntry.isDict())
        {
            WarningInFunction
                << "Ignoring non-dictionary entry "
                << dEntry << endl;
            continue;
        }

        const dictionary& dict = dEntry.dict();

        const auto action = actionNames.get("action", dict);

        // Handle manually
        if (action == topoSetSource::INVERT)
        {
            cellsToSelect.invert(mesh.nCells());
            continue;
        }

        auto source = topoSetCellSource::New
        (
            dict.get<word>("source"),
            mesh,
            dict.optionalSubDict("sourceInfo")
        );
        source->verbose(false);

        switch (action)
        {
            case topoSetSource::NEW:  // "use"
            case topoSetSource::ADD:
            case topoSetSource::SUBTRACT:
                if (topoSetSource::NEW == action)
                {
                    // "use": only use this selection (clear + ADD)
                    // NEW is handled like ADD in applyToSet()
                    cellsToSelect.reset();
                }
                source->applyToSet(action, cellsToSelect);
                break;

            case topoSetSource::SUBSET:
            {
                cellBitSet other(mesh, false);
                source->applyToSet(topoSetSource::NEW, other);

                cellsToSelect.subset(other);
            }
            break;

            default:
                // Should already have been caught
                WarningInFunction
                    << "Ignoring unhandled action '"
                    << actionNames[action] << "'" << endl;
                break;
        }
    }

    subsetter.setCellSubset(cellsToSelect.addressing());

    return true;
}


bool Foam::functionObjects::ensightWrite::update()
{
    if (meshState_ == polyMesh::UNCHANGED)
    {
        return false;
    }

    // This is heavy-handed, but with a bounding-box limited sub-mesh,
    // we don't readily know if the updates affect the subsetted mesh.

    // if (meshSubset_.hasSubMesh())
    // {
    //     ensMesh_.clear();
    //     meshSubset_.clear();
    // }
    // else if (ensMesh_)
    // {
    //     ensMesh_->expire();
    // }

    meshSubset_.clear();

    updateSubset(meshSubset_);

    meshState_ = polyMesh::UNCHANGED;

    if (!ensMesh_)
    {
        ensMesh_.reset(new ensightMesh(meshSubset_.mesh(), writeOpts_));
    }
    else if (ensMesh_->needsUpdate())
    {
        ensMesh_->correct();
    }

    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ensightWrite::readSelection(const dictionary& dict)
{
    // Ensure consistency
    ensMesh_.clear();

    meshSubset_.clear();
    meshState_ = polyMesh::TOPO_CHANGE;

    selectFields_.clear();
    dict.readEntry("fields", selectFields_);
    selectFields_.uniq();

    // Actions to define selection
    selection_ = dict.subOrEmptyDict("selection");

    return true;
}


void Foam::functionObjects::ensightWrite::updateMesh(const mapPolyMesh&)
{
    meshState_ = polyMesh::TOPO_CHANGE;
}


void Foam::functionObjects::ensightWrite::movePoints(const polyMesh&)
{
    // Only move to worse states
    if (meshState_ == polyMesh::UNCHANGED)
    {
        meshState_ = polyMesh::POINTS_MOVED;
    }
}


// ************************************************************************* //

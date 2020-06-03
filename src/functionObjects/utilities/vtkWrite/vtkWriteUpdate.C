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

#include "vtkWrite.H"
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

bool Foam::functionObjects::vtkWrite::updateSubset
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


Foam::labelList Foam::functionObjects::vtkWrite::getSelectedPatches
(
    const polyBoundaryMesh& patches
) const
{
    DynamicList<label> patchIDs(patches.size());

    for (const polyPatch& pp : patches)
    {
        if (isType<emptyPolyPatch>(pp))
        {
            continue;
        }
        else if (isA<processorPolyPatch>(pp))
        {
            break; // No processor patches
        }

        if
        (
            selectPatches_.size()
          ? selectPatches_.match(pp.name())
          : true
        )
        {
            patchIDs.append(pp.index());
        }
    }

    return patchIDs.shrink();
}


bool Foam::functionObjects::vtkWrite::update()
{
    if
    (
        meshState_ == polyMesh::UNCHANGED
     && (meshes_.size() == meshSubsets_.size())
     && (meshes_.size() == vtuMappings_.size())
    )
    {
        return false;
    }

    meshSubsets_.resize(meshes_.size());
    vtuMappings_.resize(meshes_.size());

    label regioni = 0;
    for (const word& regionName : meshes_.sortedToc())
    {
        const fvMesh& mesh = *(meshes_[regionName]);

        if (meshSubsets_.set(regioni))
        {
            meshSubsets_[regioni].clear();
        }
        else
        {
            // Mesh subsetting, or pass through
            meshSubsets_.set(regioni, new fvMeshSubset(mesh));
        }

        if (vtuMappings_.set(regioni))
        {
            // Trigger change for vtk cells too
            vtuMappings_[regioni].clear();
        }
        else
        {
            // VTU sizing and decomposition information
            vtuMappings_.set
            (
                regioni,
                new vtk::vtuCells(writeOpts_, decompose_)
            );
        }

        ++regioni;
    }

    regioni = 0;
    for (auto& subsetter : meshSubsets_)
    {
        updateSubset(subsetter);
        vtuMappings_[regioni].reset(subsetter.mesh());
        ++regioni;
    }

    meshState_ = polyMesh::UNCHANGED;
    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::vtkWrite::readSelection(const dictionary& dict)
{
    meshSubsets_.clear();
    vtuMappings_.clear();
    meshState_ = polyMesh::TOPO_CHANGE;

    // All possible meshes
    meshes_ = time_.lookupClass<fvMesh>();

    selectRegions_.clear();
    dict.readIfPresent("regions", selectRegions_);

    if (selectRegions_.empty())
    {
        selectRegions_.resize(1);
        selectRegions_.first() =
            dict.getOrDefault<word>("region", polyMesh::defaultRegion);
    }

    // Restrict to specified meshes
    meshes_.filterKeys(selectRegions_);

    if (meshes_.empty())
    {
        WarningInFunction
            << "No mesh regions selected for function object " << name()
            << nl;
    }

    selectPatches_.clear();
    dict.readIfPresent("patches", selectPatches_);

    selectFields_.clear();
    dict.readEntry("fields", selectFields_);
    selectFields_.uniq();

    // Actions to define selection
    selection_ = dict.subOrEmptyDict("selection");

    return true;
}


void Foam::functionObjects::vtkWrite::updateMesh(const mapPolyMesh&)
{
    meshState_ = polyMesh::TOPO_CHANGE;
}


void Foam::functionObjects::vtkWrite::movePoints(const polyMesh&)
{
    // Only move to worse states
    if (meshState_ == polyMesh::UNCHANGED)
    {
        meshState_ = polyMesh::POINTS_MOVED;
    }
}


// ************************************************************************* //

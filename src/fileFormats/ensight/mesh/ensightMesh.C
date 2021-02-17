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

#include "ensightMesh.H"
#include "ensightGeoFile.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::ensightMesh::internalZone = -1;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Patch names without processor patches
static wordList nonProcessorPatchNames(const polyBoundaryMesh& bmesh)
{
    #ifdef FULLDEBUG
    // Patches are output. Check that they are synced.
    bmesh.checkParallelSync(true);
    #endif

    wordList patchNames(bmesh.names());
    patchNames.resize(bmesh.nNonProcessor());

    return patchNames;
}


} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightMesh::clear()
{
    cellZoneParts_.clear();
    faceZoneParts_.clear();
    boundaryParts_.clear();
}


void Foam::ensightMesh::renumber()
{
    label partNo = 0;

    for (const label id : cellZoneParts_.sortedToc())
    {
        cellZoneParts_[id].index() = partNo++;
    }

    for (const label id : boundaryParts_.sortedToc())
    {
        boundaryParts_[id].index() = partNo++;
    }

    for (const label id : faceZoneParts_.sortedToc())
    {
        faceZoneParts_[id].index() = partNo++;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightMesh::ensightMesh
(
    const polyMesh& mesh
)
:
    ensightMesh(mesh, ensightMesh::options())
{}


Foam::ensightMesh::ensightMesh
(
    const polyMesh& mesh,
    const ensightMesh::options& opts
)
:
    options_(new options(opts)),
    mesh_(mesh),
    needsUpdate_(true)
{
    if (!option().lazy())
    {
        correct();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightMesh::correct()
{
    clear();

    const wordRes& czMatcher = option().cellZoneSelection();
    const wordRes& fzMatcher = option().faceZoneSelection();

    // Possible cellZones
    const wordList czNames =
    (
        (
            option().useCellZones()
         && (!czMatcher.empty() || option().useInternalMesh())
        )
      ? mesh_.cellZones().names()
      : wordList()
    );

    const labelList czoneIds =
    (
        czMatcher.empty()
      ? identity(czNames.size())        // All
      : czMatcher.matching(czNames)     // Selected names
    );


    // Possible faceZones
    const wordList fzNames =
    (
        option().useFaceZones()
      ? mesh_.faceZones().names()
      : wordList()
    );

    const labelList fzoneIds =
    (
        fzMatcher.empty()
      ? identity(fzNames.size())        // All
      : fzMatcher.matching(fzNames)     // Selected names
    );


    // Possible patchNames
    const wordList patchNames =
    (
        option().useBoundaryMesh()
      ? nonProcessorPatchNames(mesh_.boundaryMesh())
      : wordList()
    );

    const labelList patchIds =
    (
        option().useBoundaryMesh()
      ? stringListOps::findMatching
        (
            patchNames,
            option().patchSelection(),
            option().patchExclude()
        )
      : labelList()
    );


    // Track which cells are in a zone or not
    bitSet cellSelection;

    // Faces to be excluded from export
    bitSet excludeFace;


    // cellZones first
    for (const label zoneId : czoneIds)
    {
        const word& zoneName = czNames[zoneId];
        const cellZone& zn = mesh_.cellZones()[zoneId];

        if (returnReduce(!zn.empty(), orOp<bool>()))
        {
            // Ensure full mesh coverage
            cellSelection.resize(mesh_.nCells());

            cellSelection.set(zn);

            ensightCells& part = cellZoneParts_(zoneId);

            part.clear();
            part.identifier() = zoneId;
            part.rename(zoneName);

            part.classify(mesh_, zn);

            // Finalize
            part.reduce();
        }
    }

    if (option().useInternalMesh() && czMatcher.empty())
    {
        // The internal mesh:
        // Either the entire mesh (if no zones specified)
        // or whatever is leftover as unzoned

        if (cellZoneParts_.empty())
        {
            ensightCells& part = cellZoneParts_(internalZone);

            part.clear();
            part.identifier() = internalZone;
            part.rename("internalMesh");

            part.classify(mesh_);

            // Finalize
            part.reduce();
        }
        else
        {
            // Unzoned cells - flip selection from zoned to unzoned
            cellSelection.flip();

            if (returnReduce(cellSelection.any(), orOp<bool>()))
            {
                ensightCells& part = cellZoneParts_(internalZone);

                part.clear();
                part.identifier() = internalZone;
                part.rename("internalMesh");

                part.classify(mesh_, cellSelection);

                // Finalize
                part.reduce();
            }
        }

        // Handled all cells
        cellSelection.clearStorage();
    }
    else if (cellSelection.none())
    {
        // No internal, no cellZones selected - just ignore
        cellSelection.clearStorage();
    }


    // Face exclusion based on cellZones

    if (returnReduce(!cellSelection.empty(), orOp<bool>()))
    {
        // Ensure full mesh coverage
        excludeFace.resize(mesh_.nFaces());

        const labelList& owner = mesh_.faceOwner();

        forAll(owner, facei)
        {
            const label celli = owner[facei];

            if (!cellSelection.test(celli))
            {
                excludeFace.set(facei);
            }
        }
    }


    // Face exclusion for empty/processor types
    // so they are ignored for face zones

    if (fzoneIds.size())
    {
        // Ensure full mesh coverage
        excludeFace.resize(mesh_.nFaces());

        for (const polyPatch& p : mesh_.boundaryMesh())
        {
            const auto* procPatch = isA<processorPolyPatch>(p);

            if (isA<emptyPolyPatch>(p))
            {
                excludeFace.set(p.range());
            }
            else if (procPatch && !procPatch->owner())
            {
                // Exclude neighbour-side, retain owner-side only
                excludeFace.set(p.range());
            }
        }
    }


    // Patches
    for (const label patchId : patchIds)
    {
        const word& patchName = patchNames[patchId];
        const polyPatch& p = mesh_.boundaryMesh()[patchId];

        if (isA<emptyPolyPatch>(p))
        {
            // Skip empty patch types
            continue;
        }
        else if (isA<processorPolyPatch>(p))
        {
            // Only real (non-processor) boundaries.
            break;
        }

        ensightFaces& part = boundaryParts_(patchId);

        part.clear();
        part.identifier() = patchId;
        part.rename(patchName);

        part.classify
        (
            mesh_.faces(),
            identity(p.range()),
            boolList(),  // no flip map
            excludeFace
        );

        // Finalize
        part.reduce();

        if (!part.total())
        {
            boundaryParts_.erase(patchId);
        }
    }


    // Face zones

    for (const label zoneId : fzoneIds)
    {
        const word& zoneName = fzNames[zoneId];
        const faceZone& zn = mesh_.faceZones()[zoneId];

        ensightFaces& part = faceZoneParts_(zoneId);

        part.clear();
        part.identifier() = zoneId;
        part.rename(zoneName);

        if (zn.size())
        {
            part.classify
            (
                mesh_.faces(),
                zn,
                zn.flipMap(),
                excludeFace
            );
        }

        // Finalize
        part.reduce();

        if (!part.total())
        {
            faceZoneParts_.erase(zoneId);
        }
    }

    renumber();

    needsUpdate_ = false;
}


void Foam::ensightMesh::write
(
    ensightGeoFile& os,
    bool parallel
) const
{
    // The internalMesh, cellZones
    for (const label id : cellZoneParts_.sortedToc())
    {
        cellZoneParts_[id].write(os, mesh_, parallel);
    }

    // Patches - sorted by index
    for (const label id : boundaryParts_.sortedToc())
    {
        boundaryParts_[id].write(os, mesh_, parallel);
    }

    // Requested faceZones - sorted by index
    for (const label id : faceZoneParts_.sortedToc())
    {
        faceZoneParts_[id].write(os, mesh_, parallel);
    }
}


// ************************************************************************* //

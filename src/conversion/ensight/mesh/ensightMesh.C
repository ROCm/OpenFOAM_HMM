/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
#include "fvMesh.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "processorPolyPatch.H"
#include "mapDistribute.H"
#include "stringListOps.H"

#include "ensightFile.H"
#include "ensightGeoFile.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightMesh::clear()
{
    meshCells_.clear();
    boundaryPatchFaces_.clear();
    faceZoneFaces_.clear();
    patchLookup_.clear();
    globalPointsPtr_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightMesh::ensightMesh
(
    const fvMesh& mesh,
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


Foam::ensightMesh::ensightMesh
(
    const fvMesh& mesh,
    const IOstream::streamFormat format
)
:
    options_(new options(format)),
    mesh_(mesh),
    needsUpdate_(true)
{
    if (!option().lazy())
    {
        correct();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightMesh::~ensightMesh()
{
    deleteDemandDrivenData(options_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ensightMesh::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::ensightMesh::expire()
{
    clear();

    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


void Foam::ensightMesh::correct()
{
    clear();

    // First see if patches are allowed/disallowed
    // and if only particular patches should be selected

    label nParts = 1;  // provisionally (for internalMesh)

    if (option().usePatches())
    {
        // Patches are output. Check that they are synced.
        mesh_.boundaryMesh().checkParallelSync(true);

        wordList patchNames = mesh_.boundaryMesh().names();
        if (Pstream::parRun())
        {
            patchNames.setSize
            (
                mesh_.boundary().size()
              - mesh_.globalData().processorPatches().size()
            );
        }

        labelList matched;

        bool useAll = true;
        const wordReList& matcher = option().patchSelection();
        if (notNull(matcher))
        {
            nParts = 0; // no internalMesh

            if (!matcher.empty())
            {
                useAll = false;
                matched = findMatchingStrings
                (
                    wordRes(matcher),
                    patchNames
                );
            }
        }

        if (useAll)
        {
            matched = identity(patchNames.size());
        }

        forAll(matched, matchi)
        {
            const label patchId   = matched[matchi];
            const word& patchName = patchNames[patchId];

            // use fvPatch (not polyPatch) to automatically remove empty patches
            const fvPatch& p = mesh_.boundary()[patchId];

            // yes we most likely want this patch.
            // - can use insert or set, since hash was cleared before
            boundaryPatchFaces_.set(patchName, ensightFaces());
            ensightFaces& ensFaces = boundaryPatchFaces_[patchName];

            if (p.size())
            {
                // use local face addressing (offset = 0),
                // since this is what we'll need later when writing fields
                ensFaces.classify(p.patch());
            }
            else
            {
                // patch is empty (on this processor)
                // or the patch is 'empty' (as fvPatch type)
                ensFaces.clear();
            }

            // finalize
            ensFaces.reduce();

            if (ensFaces.total())
            {
                patchLookup_.set(patchId, patchName);
                ensFaces.index() = nParts++;
            }
            else
            {
                boundaryPatchFaces_.erase(patchName);
            }
        }

        // At this point,
        // * patchLookup_        is a map of (patchId, name)
        // * boundaryPatchFaces_ is a lookup by name for the faces elements
    }

    if (useInternalMesh())
    {
        meshCells_.index() = 0;
        meshCells_.classify(mesh_);

        // Determine parallel shared points
        globalPointsPtr_ = mesh_.globalData().mergePoints
        (
            pointToGlobal_,
            uniquePointMap_
        );
    }

    meshCells_.reduce();

    // faceZones
    if (option().useFaceZones())
    {
        // Mark boundary faces to be excluded from export
        PackedBoolList excludeFace(mesh_.nFaces()); // all false

        forAll(mesh_.boundaryMesh(), patchi)
        {
            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            if
            (
                isA<processorPolyPatch>(pp)
             && !refCast<const processorPolyPatch>(pp).owner()
            )
            {
                label bFaceI = pp.start();
                forAll(pp, i)
                {
                    excludeFace.set(bFaceI++);
                }
            }
        }

        const wordReList& matcher = option().faceZoneSelection();

        wordList selectZones = mesh_.faceZones().names();
        inplaceSubsetMatchingStrings
        (
            wordRes(matcher),
            selectZones
        );

        // have same order as later with sortedToc()
        Foam::sort(selectZones);

        // Count face types in each selected faceZone
        forAll(selectZones, zonei)
        {
            const word& zoneName = selectZones[zonei];
            const label zoneID = mesh_.faceZones().findZoneID(zoneName);
            const faceZone& fz = mesh_.faceZones()[zoneID];

            // yes we most likely want this zone
            // - can use insert or set, since hash was cleared before
            faceZoneFaces_.set(zoneName, ensightFaces());
            ensightFaces& ensFaces = faceZoneFaces_[zoneName];

            if (fz.size())
            {
                ensFaces.classify
                (
                    mesh_.faces(),
                    fz,
                    fz.flipMap(),
                    excludeFace
                );
            }

            // finalize
            ensFaces.reduce();

            if (ensFaces.total())
            {
                ensFaces.index() = nParts++;
            }
            else
            {
                faceZoneFaces_.erase(zoneName);
            }
        }
    }

    needsUpdate_ = false;
}


void Foam::ensightMesh::write(ensightGeoFile& os) const
{
    if (useInternalMesh())
    {
        label nPoints = globalPoints().size();

        const pointField uniquePoints(mesh_.points(), uniquePointMap_);

        // writePartHeader(os, 0, "internalMesh");
        // beginCoordinates(os, nPoints);
        writeAllPoints
        (
            meshCells_.index(),
            "internalMesh",
            nPoints,
            uniquePoints,
            os
        );

        writeCellConnectivity(meshCells_, pointToGlobal_, os);
    }


    //
    // write patches
    // use sortedToc for extra safety
    //
    const labelList patchIds = patchLookup_.sortedToc();
    forAll(patchIds, listi)
    {
        const label patchId   = patchIds[listi];
        const word& patchName = patchLookup_[patchId];
        const ensightFaces& ensFaces = boundaryPatchFaces_[patchName];

        const polyPatch& pp = mesh_.boundaryMesh()[patchId];

        // Renumber the patch points/faces into unique points
        labelList pointToGlobal;
        labelList uniqueMeshPointLabels;
        autoPtr<globalIndex> globalPointsPtr =
            mesh_.globalData().mergePoints
            (
                pp.meshPoints(),
                pp.meshPointMap(),
                pointToGlobal, // local point to unique global index
                uniqueMeshPointLabels // unique global points
            );

        // Renumber the patch faces,
        // from local patch indexing to unique global index
        faceList patchFaces(pp.localFaces());
        forAll(patchFaces, i)
        {
            inplaceRenumber(pointToGlobal, patchFaces[i]);
        }

        writeAllPoints
        (
            ensFaces.index(),
            patchName,
            globalPointsPtr().size(),
            pointField(mesh_.points(), uniqueMeshPointLabels),
            os
        );

        writeFaceConnectivity(ensFaces, patchFaces, os);
    }


    //
    // write faceZones, if requested
    //
    const wordList zoneNames = faceZoneFaces_.sortedToc();
    forAll(zoneNames, zonei)
    {
        const word& zoneName = zoneNames[zonei];
        const ensightFaces& ensFaces = faceZoneFaces_[zoneName];

        // Use the properly sorted faceIds (ensightFaces) and do NOT use the
        // faceZone directly, otherwise the point-maps will not correspond.
        // - perform face-flipping later

        indirectPrimitivePatch pp
        (
            IndirectList<face>(mesh_.faces(), ensFaces.faceIds()),
            mesh_.points()
        );

        // Renumber the points/faces into unique points
        labelList pointToGlobal;
        labelList uniqueMeshPointLabels;
        autoPtr<globalIndex> globalPointsPtr =
            mesh_.globalData().mergePoints
            (
                pp.meshPoints(),
                pp.meshPointMap(),
                pointToGlobal, // local point to unique global index
                uniqueMeshPointLabels // unique global points
            );

        // Renumber the faces belonging to the faceZone,
        // from local numbering to unique global index.
        // Also a good place to perform face flipping
        const boolList& flip = ensFaces.flipMap();
        faceList patchFaces(pp.localFaces());
        forAll(patchFaces, facei)
        {
            face& f = patchFaces[facei];

            if (flip[facei])
            {
                f.flip();
            }

            inplaceRenumber(pointToGlobal, f);
        }

        writeAllPoints
        (
            ensFaces.index(),
            zoneName,
            globalPointsPtr().size(),
            pointField(mesh_.points(), uniqueMeshPointLabels),
            os
        );

        writeFaceConnectivity(ensFaces, patchFaces, os, true);
    }
}


// ************************************************************************* //

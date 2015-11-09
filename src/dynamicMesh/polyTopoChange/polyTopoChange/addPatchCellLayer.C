/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "addPatchCellLayer.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "syncTools.H"
#include "polyAddPoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyAddCell.H"
#include "globalIndex.H"
#include "PatchTools.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(addPatchCellLayer, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::addPatchCellLayer::nbrFace
(
    const labelListList& edgeFaces,
    const label edgeI,
    const label faceI
)
{
    const labelList& eFaces = edgeFaces[edgeI];

    if (eFaces.size() == 2)
    {
        return (eFaces[0] != faceI ? eFaces[0] : eFaces[1]);
    }
    else
    {
        return -1;
    }
}


void Foam::addPatchCellLayer::addVertex
(
    const label pointI,
    face& f,
    label& fp
)
{
    if (fp == 0)
    {
        f[fp++] = pointI;
    }
    else
    {
        if (f[fp-1] != pointI && f[0] != pointI)
        {
            f[fp++] = pointI;
        }
    }
}


// Is edge to the same neighbour? (and needs extrusion and has not been
// dealt with already)
bool Foam::addPatchCellLayer::sameEdgeNeighbour
(
    const indirectPrimitivePatch& pp,
    const labelListList& globalEdgeFaces,
    const boolList& doneEdge,
    const label thisGlobalFaceI,
    const label nbrGlobalFaceI,
    const label edgeI
) const
{
    const edge& e = pp.edges()[edgeI];

    return
        !doneEdge[edgeI]                            // not yet handled
     && (
            addedPoints_[e[0]].size()               // is extruded
         || addedPoints_[e[1]].size()
        )
     && (
            nbrFace(globalEdgeFaces, edgeI, thisGlobalFaceI)
         == nbrGlobalFaceI  // is to same neighbour
        );
}


// Collect consecutive string of edges that connects the same two
// (possibly coupled) faces. Returns -1 if no unvisited edge can be found.
// Otherwise returns start and end index in face.
Foam::labelPair Foam::addPatchCellLayer::getEdgeString
(
    const indirectPrimitivePatch& pp,
    const labelListList& globalEdgeFaces,
    const boolList& doneEdge,
    const label patchFaceI,
    const label globalFaceI
) const
{
    const labelList& fEdges = pp.faceEdges()[patchFaceI];

    label startFp = -1;
    label endFp = -1;

    // Get edge that hasn't been done yet but needs extrusion
    forAll(fEdges, fp)
    {
        label edgeI = fEdges[fp];
        const edge& e = pp.edges()[edgeI];

        if
        (
            !doneEdge[edgeI]
         && ( addedPoints_[e[0]].size() || addedPoints_[e[1]].size() )
        )
        {
            startFp = fp;
            break;
        }
    }

    if (startFp != -1)
    {
        // We found an edge that needs extruding but hasn't been done yet.
        // Now find the face on the other side
        label nbrGlobalFaceI = nbrFace
        (
            globalEdgeFaces,
            fEdges[startFp],
            globalFaceI
        );

        if (nbrGlobalFaceI == -1)
        {
            // Proper boundary edge. Only extrude single edge.
            endFp = startFp;
        }
        else
        {
            // Search back for edge
            // - which hasn't been handled yet
            // - with same neighbour
            // - that needs extrusion
            while (true)
            {
                label prevFp = fEdges.rcIndex(startFp);

                if
                (
                    !sameEdgeNeighbour
                    (
                        pp,
                        globalEdgeFaces,
                        doneEdge,
                        globalFaceI,
                        nbrGlobalFaceI,
                        fEdges[prevFp]
                    )
                )
                {
                    break;
                }
                startFp = prevFp;
            }

            // Search forward for end of string
            endFp = startFp;
            while (true)
            {
                label nextFp = fEdges.fcIndex(endFp);

                if
                (
                    !sameEdgeNeighbour
                    (
                        pp,
                        globalEdgeFaces,
                        doneEdge,
                        globalFaceI,
                        nbrGlobalFaceI,
                        fEdges[nextFp]
                    )
                )
                {
                    break;
                }
                endFp = nextFp;
            }
        }
    }

    return labelPair(startFp, endFp);
}


Foam::label Foam::addPatchCellLayer::addSideFace
(
    const indirectPrimitivePatch& pp,
    const labelListList& addedCells,    // per pp face the new extruded cell
    const face& newFace,
    const label newPatchID,
    const label zoneI,
    const bool edgeFlip,
    const label inflateFaceI,

    const label ownFaceI,               // pp face that provides owner
    const label nbrFaceI,
    const label meshEdgeI,              // corresponding mesh edge
    const label layerI,                 // layer
    const label numEdgeFaces,           // number of layers for edge
    const labelList& meshFaces,         // precalculated edgeFaces
    polyTopoChange& meshMod
) const
{
    label addedFaceI = -1;


    // Is patch edge external edge of indirectPrimitivePatch?
    if (nbrFaceI == -1)
    {
        // External edge so external face.

        // Determine if different number of layer on owner and neighbour side
        // (relevant only for coupled faces). See section for internal edge
        // below.

        label layerOwn;

        if (addedCells[ownFaceI].size() < numEdgeFaces)
        {
            label offset = numEdgeFaces - addedCells[ownFaceI].size();
            if (layerI <= offset)
            {
                layerOwn = 0;
            }
            else
            {
                layerOwn = layerI - offset;
            }
        }
        else
        {
            layerOwn = layerI;
        }


        //Pout<< "Added boundary face:" << newFace
        //    << " own:" << addedCells[ownFaceI][layerOwn]
        //    << " patch:" << newPatchID
        //    << endl;

        addedFaceI = meshMod.setAction
        (
            polyAddFace
            (
                newFace,                        // face
                addedCells[ownFaceI][layerOwn], // owner
                -1,                             // neighbour
                -1,                             // master point
                -1,                             // master edge
                inflateFaceI,                   // master face
                false,                          // flux flip
                newPatchID,                     // patch for face
                zoneI,                          // zone for face
                edgeFlip                        // face zone flip
            )
        );
    }
    else
    {
        // When adding side faces we need to modify neighbour and owners
        // in region where layer mesh is stopped. Determine which side
        // has max number of faces and make sure layers match closest to
        // original pp if there are different number of layers.

        label layerNbr;
        label layerOwn;

        if (addedCells[ownFaceI].size() > addedCells[nbrFaceI].size())
        {
            label offset =
                addedCells[ownFaceI].size() - addedCells[nbrFaceI].size();

            layerOwn = layerI;

            if (layerI <= offset)
            {
                layerNbr = 0;
            }
            else
            {
                layerNbr = layerI - offset;
            }
        }
        else if (addedCells[nbrFaceI].size() > addedCells[ownFaceI].size())
        {
            label offset =
                addedCells[nbrFaceI].size() - addedCells[ownFaceI].size();

            layerNbr = layerI;

            if (layerI <= offset)
            {
                layerOwn = 0;
            }
            else
            {
                layerOwn = layerI - offset;
            }
        }
        else
        {
            // Same number of layers on both sides.
            layerNbr = layerI;
            layerOwn = layerI;
        }


        // Check mesh internal faces using edge to initialise
        label inflateEdgeI = -1;
        if (addToMesh_)
        {
            forAll(meshFaces, i)
            {
                if (mesh_.isInternalFace(meshFaces[i]))
                {
                    // meshEdge uses internal faces so ok to inflate from it
                    inflateEdgeI = meshEdgeI;
                    break;
                }
            }
        }


        addedFaceI = meshMod.setAction
        (
            polyAddFace
            (
                newFace,                        // face
                addedCells[ownFaceI][layerOwn], // owner
                addedCells[nbrFaceI][layerNbr], // neighbour
                -1,                             // master point
                inflateEdgeI,                   // master edge
                -1,                             // master face
                false,                          // flux flip
                -1,                             // patch for face
                zoneI,                          // zone for face
                edgeFlip                        // face zone flip
            )
        );

       //Pout<< "Added internal face:" << newFace
        //    << " own:" << addedCells[ownFaceI][layerOwn]
        //    << " nei:" << addedCells[nbrFaceI][layerNbr]
        //    << endl;
    }

    return addedFaceI;
}


Foam::label Foam::addPatchCellLayer::findProcPatch
(
    const polyMesh& mesh,
    const label nbrProcID
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(mesh.globalData().processorPatches(), i)
    {
        label patchI = mesh.globalData().processorPatches()[i];

        if
        (
            refCast<const processorPolyPatch>(patches[patchI]).neighbProcNo()
         == nbrProcID
        )
        {
            return patchI;
        }
    }
    return -1;
}


void Foam::addPatchCellLayer::setFaceProps
(
    const polyMesh& mesh,
    const label faceI,

    label& patchI,
    label& zoneI,
    bool& zoneFlip
)
{
    patchI = mesh.boundaryMesh().whichPatch(faceI);
    zoneI = mesh.faceZones().whichZone(faceI);
    if (zoneI != -1)
    {
        label index = mesh.faceZones()[zoneI].whichFace(faceI);
        zoneFlip = mesh.faceZones()[zoneI].flipMap()[index];
    }
}


void Foam::addPatchCellLayer::setFaceProps
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    const label ppEdgeI,
    const label faceI,

    label& patchI,
    label& zoneI,
    bool& zoneFlip,
    label& inflateFaceI
)
{
    setFaceProps
    (
        mesh,
        faceI,

        patchI,
        zoneI,
        zoneFlip
    );

    if (patchI != -1 || zoneI != -1)
    {
        inflateFaceI = faceI;
    }

    if (zoneI != -1)
    {
        // Correct flip for patch edge ordering
        const edge& ppEdge = pp.edges()[ppEdgeI];
        const edge patchEdge
        (
            pp.meshPoints()[ppEdge[0]],
            pp.meshPoints()[ppEdge[1]]
        );

        const face& f = mesh.faces()[faceI];
        bool found = false;
        forAll(f, fp)
        {
            const edge e(f[fp], f.nextLabel(fp));
            int stat = edge::compare(e, patchEdge);
            if (stat == 1)
            {
                found = true;
                break;
            }
            else if (stat == -1)
            {
                found = true;
                zoneFlip = !zoneFlip;
                break;
            }
        }

        if (!found)
        {
            FatalErrorIn("addPatchCellLayer::setFaceProps(..)")
                << "Problem: cannot find patch edge " << ppEdgeI
                << " with mesh vertices " << patchEdge
                << " at " << patchEdge.line(mesh.points())
                << " is not in face " << faceI << " with mesh vertices "
                << f
                << exit(FatalError);
        }
    }
}


void Foam::addPatchCellLayer::findZoneFace
(
    const bool useInternalFaces,
    const bool useBoundaryFaces,

    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    const label ppEdgeI,
    const UIndirectList<label>& excludeFaces,
    const labelList& meshFaces,

    label& inflateFaceI,
    label& patchI,
    label& zoneI,
    bool& zoneFlip
)
{
    inflateFaceI = -1;
    patchI = -1;
    zoneI = -1;
    zoneFlip = false;

    forAll(meshFaces, k)
    {
        label faceI = meshFaces[k];

        if
        (
            (findIndex(excludeFaces, faceI) == -1)
         && (
                (mesh.isInternalFace(faceI) && useInternalFaces)
             || (!mesh.isInternalFace(faceI) && useBoundaryFaces)
            )
        )
        {
            setFaceProps
            (
                mesh,
                pp,
                ppEdgeI,
                faceI,

                patchI,
                zoneI,
                zoneFlip,
                inflateFaceI
            );

            if (zoneI != -1 || patchI != -1)
            {
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::addPatchCellLayer::addPatchCellLayer
(
    const polyMesh& mesh,
    const bool addToMesh
)
:
    mesh_(mesh),
    addToMesh_(addToMesh),
    addedPoints_(0),
    layerFaces_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelListList Foam::addPatchCellLayer::addedCells
(
    const polyMesh& mesh,
    const labelListList& layerFaces
)
{
    labelListList layerCells(layerFaces.size());

    forAll(layerFaces, patchFaceI)
    {
        const labelList& faceLabels = layerFaces[patchFaceI];

        if (faceLabels.size())
        {
            labelList& added = layerCells[patchFaceI];
            added.setSize(faceLabels.size()-1);

            for (label i = 0; i < faceLabels.size()-1; i++)
            {
                added[i] = mesh.faceNeighbour()[faceLabels[i]];
            }
        }
    }
    return layerCells;
}


Foam::labelListList Foam::addPatchCellLayer::addedCells() const
{
    return addedCells(mesh_, layerFaces_);
}


// Calculate global faces per pp edge.
Foam::labelListList Foam::addPatchCellLayer::globalEdgeFaces
(
    const polyMesh& mesh,
    const globalIndex& globalFaces,
    const indirectPrimitivePatch& pp
)
{
    // Precalculate mesh edges for pp.edges.
    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    // From mesh edge to global face labels. Non-empty sublists only for
    // pp edges.
    labelListList globalEdgeFaces(mesh.nEdges());

    const labelListList& edgeFaces = pp.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];

        const labelList& eFaces = edgeFaces[edgeI];

        // Store face and processor as unique tag.
        labelList& globalEFaces = globalEdgeFaces[meshEdgeI];
        globalEFaces.setSize(eFaces.size());
        forAll(eFaces, i)
        {
            globalEFaces[i] = globalFaces.toGlobal(pp.addressing()[eFaces[i]]);
        }
    }

    // Synchronise across coupled edges.
    syncTools::syncEdgeList
    (
        mesh,
        globalEdgeFaces,
        uniqueEqOp(),
        labelList()             // null value
    );

    // Extract pp part
    return labelListList(UIndirectList<labelList>(globalEdgeFaces, meshEdges));
}


void Foam::addPatchCellLayer::calcExtrudeInfo
(
    const bool zoneFromAnyFace,

    const polyMesh& mesh,
    const globalIndex& globalFaces,
    const labelListList& globalEdgeFaces,
    const indirectPrimitivePatch& pp,

    labelList& edgePatchID,
    label& nPatches,
    Map<label>& nbrProcToPatch,
    Map<label>& patchToNbrProc,
    labelList& edgeZoneID,
    boolList& edgeFlip,
    labelList& inflateFaceID
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const globalMeshData& gd = mesh.globalData();

    // Precalculate mesh edges for pp.edges.
    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    edgePatchID.setSize(pp.nEdges());
    edgePatchID = -1;
    edgeZoneID.setSize(pp.nEdges());
    edgeZoneID = -1;
    edgeFlip.setSize(pp.nEdges());
    edgeFlip = false;


    // Determine properties for faces extruded from edges
    // - edge inbetween two different processors:
    //      - extrude as patch face on correct processor
    // - edge at end of patch (so edgeFaces.size() == 1):
    //      - use mesh face that is a boundary face
    // - edge internal to patch (so edgeFaces.size() == 2):



    // These also get determined but not (yet) exported:
    // - whether face is created from other face or edge

    inflateFaceID.setSize(pp.nEdges(), -1);

    nPatches = patches.size();



    // Pass1:
    // For all edges inbetween two processors: see if matches to existing
    // processor patch or create interprocessor-patch if necessary.
    // Sets edgePatchID[edgeI] but none of the other quantities

    forAll(globalEdgeFaces, edgeI)
    {
        const labelList& eGlobalFaces = globalEdgeFaces[edgeI];
        if
        (
            eGlobalFaces.size() == 2
         && pp.edgeFaces()[edgeI].size() == 1
        )
        {
            // Locally but not globally a boundary edge. Hence a coupled
            // edge. Find the patch to use if on different processors.

            label f0 = eGlobalFaces[0];
            label f1 = eGlobalFaces[1];

            label otherProcI = -1;
            if (globalFaces.isLocal(f0) && !globalFaces.isLocal(f1))
            {
                otherProcI = globalFaces.whichProcID(f1);
            }
            else if (!globalFaces.isLocal(f0) && globalFaces.isLocal(f1))
            {
                otherProcI = globalFaces.whichProcID(f0);
            }


            if (otherProcI != -1)
            {
                if (findIndex(gd[Pstream::myProcNo()], otherProcI) != -1)
                {
                    // There is already a processorPolyPatch to otherProcI.
                    // Use it. Note that we can only index procPatchMap
                    // if the processor actually is a neighbour processor
                    // so that is why we first check.
                    edgePatchID[edgeI] = gd.procPatchMap()[otherProcI];
                }
                else
                {
                    // Cannot find a patch to processor. See if already
                    // marked for addition
                    if (nbrProcToPatch.found(otherProcI))
                    {
                        edgePatchID[edgeI] = nbrProcToPatch[otherProcI];
                    }
                    else
                    {
                        edgePatchID[edgeI] = nPatches;
                        nbrProcToPatch.insert(otherProcI, nPatches);
                        patchToNbrProc.insert(nPatches, otherProcI);
                        nPatches++;
                    }
                }
            }
        }
    }


    // Pass2: determine face properties for all other edges
    // ----------------------------------------------------

    const labelListList& edgeFaces = pp.edgeFaces();

    DynamicList<label> dynMeshEdgeFaces;

    forAll(edgeFaces, edgeI)
    {
        if (edgePatchID[edgeI] == -1)
        {
            UIndirectList<label> ppFaces(pp.addressing(), edgeFaces[edgeI]);

            label meshEdgeI = meshEdges[edgeI];
            const labelList& meshFaces = mesh.edgeFaces
            (
                meshEdgeI,
                dynMeshEdgeFaces
            );

            if (edgeFaces[edgeI].size() == 2)
            {
                // Internal edge. Look at any face (internal or boundary) to
                // determine extrusion properties. First one that has zone
                // info wins

                label dummyPatchI = -1;
                findZoneFace
                (
                    true,               // useInternalFaces,
                    zoneFromAnyFace,    // useBoundaryFaces,

                    mesh,
                    pp,
                    edgeI,


                    ppFaces,            // excludeFaces,
                    meshFaces,          // meshFaces,

                    inflateFaceID[edgeI],
                    dummyPatchI,        // do not use patch info
                    edgeZoneID[edgeI],
                    edgeFlip[edgeI]
                );
            }
            else
            {
                // Proper, uncoupled patch edge

                findZoneFace
                (
                    false,      // useInternalFaces,
                    true,       // useBoundaryFaces,

                    mesh,
                    pp,
                    edgeI,


                    ppFaces,    // excludeFaces,
                    meshFaces,  // meshFaces,

                    inflateFaceID[edgeI],
                    edgePatchID[edgeI],
                    edgeZoneID[edgeI],
                    edgeFlip[edgeI]
                );
            }
        }
    }



    // Now hopefully every boundary edge has a edge patch. Check
    if (debug)
    {
        forAll(edgeFaces, edgeI)
        {
            if (edgeFaces[edgeI].size() == 1 && edgePatchID[edgeI] == -1)
            {
                const edge& e = pp.edges()[edgeI];
                //FatalErrorIn("addPatchCellLayer::calcExtrudeInfo(..)")
                WarningIn("addPatchCellLayer::calcExtrudeInfo(..)")
                    << "Have no edgePatchID for edge " << edgeI << " points "
                    << pp.points()[pp.meshPoints()[e[0]]]
                    << pp.points()[pp.meshPoints()[e[1]]]
                    //<< abort(FatalError);
                    << endl;
            }
        }
    }



    // Pass3: for any faces set in pass1 see if we can find a processor face
    //        to inherit from (we only have a patch, not a patch face)
    forAll(edgeFaces, edgeI)
    {
        if
        (
            edgeFaces[edgeI].size() == 1
         && edgePatchID[edgeI] != -1
         && inflateFaceID[edgeI] == -1
        )
        {
            // 1. Do we have a boundary face to inflate from

            label myFaceI = pp.addressing()[edgeFaces[edgeI][0]];

            // Pick up any boundary face on this edge and use its properties
            label meshEdgeI = meshEdges[edgeI];
            const labelList& meshFaces = mesh.edgeFaces
            (
                meshEdgeI,
                dynMeshEdgeFaces
            );

            forAll(meshFaces, k)
            {
                label faceI = meshFaces[k];

                if (faceI != myFaceI && !mesh.isInternalFace(faceI))
                {
                    if (patches.whichPatch(faceI) == edgePatchID[edgeI])
                    {
                        setFaceProps
                        (
                            mesh,
                            pp,
                            edgeI,
                            faceI,

                            edgePatchID[edgeI],
                            edgeZoneID[edgeI],
                            edgeFlip[edgeI],
                            inflateFaceID[edgeI]
                        );
                        break;
                    }
                }
            }
        }
    }



    // Sync all data:
    // - edgePatchID : might be local in case of processor patch. Do not
    //   sync for now
    // - inflateFaceID: local. Do not sync
    // - edgeZoneID : global. sync.
    // - edgeFlip : global. sync.

    if (Pstream::parRun())
    {
        const globalMeshData& gd = mesh.globalData();
        const indirectPrimitivePatch& cpp = gd.coupledPatch();

        labelList patchEdges;
        labelList coupledEdges;
        PackedBoolList sameEdgeOrientation;
        PatchTools::matchEdges
        (
            pp,
            cpp,
            patchEdges,
            coupledEdges,
            sameEdgeOrientation
        );

        // Convert data on pp edges to data on coupled patch
        labelList cppEdgeZoneID(cpp.nEdges(), -1);
        boolList cppEdgeFlip(cpp.nEdges(), false);
        forAll(coupledEdges, i)
        {
            label cppEdgeI = coupledEdges[i];
            label ppEdgeI = patchEdges[i];

            cppEdgeZoneID[cppEdgeI] = edgeZoneID[ppEdgeI];
            if (sameEdgeOrientation[i])
            {
                cppEdgeFlip[cppEdgeI] = edgeFlip[ppEdgeI];
            }
            else
            {
                cppEdgeFlip[cppEdgeI] = !edgeFlip[ppEdgeI];
            }
        }

        // Sync
        const globalIndexAndTransform& git = gd.globalTransforms();
        const mapDistribute& edgeMap = gd.globalEdgeSlavesMap();

        globalMeshData::syncData
        (
            cppEdgeZoneID,
            gd.globalEdgeSlaves(),
            gd.globalEdgeTransformedSlaves(),
            edgeMap,
            git,
            maxEqOp<label>(),
            dummyTransform()
        );
        globalMeshData::syncData
        (
            cppEdgeFlip,
            gd.globalEdgeSlaves(),
            gd.globalEdgeTransformedSlaves(),
            edgeMap,
            git,
            andEqOp<bool>(),
            dummyTransform()
        );

        // Convert data on coupled edges to pp edges
        forAll(coupledEdges, i)
        {
            label cppEdgeI = coupledEdges[i];
            label ppEdgeI = patchEdges[i];

            edgeZoneID[ppEdgeI] = cppEdgeZoneID[cppEdgeI];
            if (sameEdgeOrientation[i])
            {
                edgeFlip[ppEdgeI] = cppEdgeFlip[cppEdgeI];
            }
            else
            {
                edgeFlip[ppEdgeI] = !cppEdgeFlip[cppEdgeI];
            }
        }
    }
}


void Foam::addPatchCellLayer::setRefinement
(
    const globalIndex& globalFaces,
    const labelListList& globalEdgeFaces,
    const scalarField& expansionRatio,
    const indirectPrimitivePatch& pp,

    const labelList& edgePatchID,
    const labelList& edgeZoneID,
    const boolList& edgeFlip,
    const labelList& inflateFaceID,

    const labelList& exposedPatchID,
    const labelList& nFaceLayers,
    const labelList& nPointLayers,
    const vectorField& firstLayerDisp,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "addPatchCellLayer::setRefinement : Adding up to "
            << gMax(nPointLayers)
            << " layers of cells to indirectPrimitivePatch with "
            << pp.nPoints() << " points" << endl;
    }

    if
    (
        pp.nPoints() != firstLayerDisp.size()
     || pp.nPoints() != nPointLayers.size()
     || pp.size() != nFaceLayers.size()
    )
    {
        FatalErrorIn
        (
            "addPatchCellLayer::setRefinement"
            "(const scalar, const indirectPrimitivePatch&"
            ", const labelList&, const vectorField&, polyTopoChange&)"
        )   << "Size of new points is not same as number of points used by"
            << " the face subset" << endl
            << "  patch.nPoints:" << pp.nPoints()
            << "  displacement:" << firstLayerDisp.size()
            << "  nPointLayers:" << nPointLayers.size() << nl
            << " patch.nFaces:" << pp.size()
            << "  nFaceLayers:" << nFaceLayers.size()
            << abort(FatalError);
    }

    forAll(nPointLayers, i)
    {
        if (nPointLayers[i] < 0)
        {
            FatalErrorIn
            (
                "addPatchCellLayer::setRefinement"
                "(const scalar, const indirectPrimitivePatch&"
                ", const labelList&, const vectorField&, polyTopoChange&)"
            )   << "Illegal number of layers " << nPointLayers[i]
                << " at patch point " << i << abort(FatalError);
        }
    }
    forAll(nFaceLayers, i)
    {
        if (nFaceLayers[i] < 0)
        {
            FatalErrorIn
            (
                "addPatchCellLayer::setRefinement"
                "(const scalar, const indirectPrimitivePatch&"
                ", const labelList&, const vectorField&, polyTopoChange&)"
            )   << "Illegal number of layers " << nFaceLayers[i]
                << " at patch face " << i << abort(FatalError);
        }
    }

    forAll(globalEdgeFaces, edgeI)
    {
        if (globalEdgeFaces[edgeI].size() > 2)
        {
            const edge& e = pp.edges()[edgeI];

            if (nPointLayers[e[0]] > 0 || nPointLayers[e[1]] > 0)
            {
                FatalErrorIn
                (
                    "addPatchCellLayer::setRefinement"
                    "(const scalar, const indirectPrimitivePatch&"
                    ", const labelList&, const vectorField&, polyTopoChange&)"
                )   << "Trying to extrude edge "
                    << e.line(pp.localPoints())
                    << " which is non-manifold (has "
                    << globalEdgeFaces[edgeI].size()
                    << " local or coupled faces using it)"
                    << " of which "
                    << pp.edgeFaces()[edgeI].size()
                    << " local"
                    << abort(FatalError);
            }
        }
    }


    const labelList& meshPoints = pp.meshPoints();

    // Some storage for edge-face-addressing.
    DynamicList<label> ef;

    // Precalculate mesh edges for pp.edges.
    const labelList meshEdges(pp.meshEdges(mesh_.edges(), mesh_.pointEdges()));

    if (debug)
    {
        // Check synchronisation
        // ~~~~~~~~~~~~~~~~~~~~~

        {
            labelList n(mesh_.nPoints(), 0);
            UIndirectList<label>(n, meshPoints) = nPointLayers;
            syncTools::syncPointList(mesh_, n, maxEqOp<label>(), label(0));

            // Non-synced
            forAll(meshPoints, i)
            {
                label meshPointI = meshPoints[i];

                if (n[meshPointI] != nPointLayers[i])
                {
                    FatalErrorIn
                    (
                        "addPatchCellLayer::setRefinement"
                        "(const scalar, const indirectPrimitivePatch&"
                        ", const labelList&, const vectorField&"
                        ", polyTopoChange&)"
                    )   << "At mesh point:" << meshPointI
                        << " coordinate:" << mesh_.points()[meshPointI]
                        << " specified nLayers:" << nPointLayers[i] << endl
                        << "On coupled point a different nLayers:"
                        << n[meshPointI] << " was specified."
                        << abort(FatalError);
                }
            }


            // Check that nPointLayers equals the max layers of connected faces
            // (or 0). Anything else makes no sense.
            labelList nFromFace(mesh_.nPoints(), 0);
            forAll(nFaceLayers, i)
            {
                const face& f = pp[i];

                forAll(f, fp)
                {
                    label pointI = f[fp];

                    nFromFace[pointI] = max(nFromFace[pointI], nFaceLayers[i]);
                }
            }
            syncTools::syncPointList
            (
                mesh_,
                nFromFace,
                maxEqOp<label>(),
                label(0)
            );

            forAll(nPointLayers, i)
            {
                label meshPointI = meshPoints[i];

                if
                (
                    nPointLayers[i] > 0
                 && nPointLayers[i] != nFromFace[meshPointI]
                )
                {
                    FatalErrorIn
                    (
                        "addPatchCellLayer::setRefinement"
                        "(const scalar, const indirectPrimitivePatch&"
                        ", const labelList&, const vectorField&"
                        ", polyTopoChange&)"
                    )   << "At mesh point:" << meshPointI
                        << " coordinate:" << mesh_.points()[meshPointI]
                        << " specified nLayers:" << nPointLayers[i] << endl
                        << "but the max nLayers of surrounding faces is:"
                        << nFromFace[meshPointI]
                        << abort(FatalError);
                }
            }
        }

        {
            pointField d(mesh_.nPoints(), vector::max);
            UIndirectList<point>(d, meshPoints) = firstLayerDisp;
            syncTools::syncPointList
            (
                mesh_,
                d,
                minEqOp<vector>(),
                vector::max
            );

            forAll(meshPoints, i)
            {
                label meshPointI = meshPoints[i];

                if (mag(d[meshPointI] - firstLayerDisp[i]) > SMALL)
                {
                    FatalErrorIn
                    (
                        "addPatchCellLayer::setRefinement"
                        "(const scalar, const indirectPrimitivePatch&"
                        ", const labelList&, const vectorField&"
                        ", polyTopoChange&)"
                    )   << "At mesh point:" << meshPointI
                        << " coordinate:" << mesh_.points()[meshPointI]
                        << " specified displacement:" << firstLayerDisp[i]
                        << endl
                        << "On coupled point a different displacement:"
                        << d[meshPointI] << " was specified."
                        << abort(FatalError);
                }
            }
        }

        // Check that edges of pp (so ones that become boundary faces)
        // connect to only one boundary face. Guarantees uniqueness of
        // patch that they go into so if this is a coupled patch both
        // sides decide the same.
        // ~~~~~~~~~~~~~~~~~~~~~~

        for (label edgeI = pp.nInternalEdges(); edgeI < pp.nEdges(); edgeI++)
        {
            const edge& e = pp.edges()[edgeI];

            if (nPointLayers[e[0]] > 0 || nPointLayers[e[1]] > 0)
            {
                // Edge is to become a face

                const labelList& eFaces = pp.edgeFaces()[edgeI];

                // First check: pp should be single connected.
                if (eFaces.size() != 1)
                {
                    FatalErrorIn
                    (
                        "addPatchCellLayer::setRefinement"
                        "(const scalar, const indirectPrimitivePatch&"
                        ", const labelList&, const vectorField&"
                        ", polyTopoChange&)"
                    )   << "boundary-edge-to-be-extruded:"
                        << pp.points()[meshPoints[e[0]]]
                        << pp.points()[meshPoints[e[1]]]
                        << " has more than two faces using it:" << eFaces
                        << abort(FatalError);
                }

                label myFaceI = pp.addressing()[eFaces[0]];

                label meshEdgeI = meshEdges[edgeI];

                // Mesh faces using edge
                const labelList& meshFaces = mesh_.edgeFaces(meshEdgeI, ef);

                // Check that there is only one patchface using edge.
                const polyBoundaryMesh& patches = mesh_.boundaryMesh();

                label bFaceI = -1;

                forAll(meshFaces, i)
                {
                    label faceI = meshFaces[i];

                    if (faceI != myFaceI)
                    {
                        if (!mesh_.isInternalFace(faceI))
                        {
                            if (bFaceI == -1)
                            {
                                bFaceI = faceI;
                            }
                            else
                            {
                                FatalErrorIn
                                (
                                    "addPatchCellLayer::setRefinement"
                                    "(const scalar"
                                    ", const indirectPrimitivePatch&"
                                    ", const labelList&, const vectorField&"
                                    ", polyTopoChange&)"
                                )   << "boundary-edge-to-be-extruded:"
                                    << pp.points()[meshPoints[e[0]]]
                                    << pp.points()[meshPoints[e[1]]]
                                    << " has more than two boundary faces"
                                    << " using it:"
                                    << bFaceI << " fc:"
                                    << mesh_.faceCentres()[bFaceI]
                                    << " patch:" << patches.whichPatch(bFaceI)
                                    << " and " << faceI << " fc:"
                                    << mesh_.faceCentres()[faceI]
                                    << " patch:" << patches.whichPatch(faceI)
                                    << abort(FatalError);
                            }
                        }
                    }
                }
            }
        }
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Precalculated patchID for each patch face
    labelList patchID(pp.size());

    forAll(pp, patchFaceI)
    {
        label meshFaceI = pp.addressing()[patchFaceI];

        patchID[patchFaceI] = patches.whichPatch(meshFaceI);
    }


    // From master point (in patch point label) to added points (in mesh point
    // label)
    addedPoints_.setSize(pp.nPoints());

    // Mark points that do not get extruded by setting size of addedPoints_ to 0
    label nTruncated = 0;

    forAll(nPointLayers, patchPointI)
    {
        if (nPointLayers[patchPointI] > 0)
        {
            addedPoints_[patchPointI].setSize(nPointLayers[patchPointI]);
        }
        else
        {
            nTruncated++;
        }
    }

    if (debug)
    {
        Pout<< "Not adding points at " << nTruncated << " out of "
            << pp.nPoints() << " points" << endl;
    }


    //
    // Create new points
    //

    // If creating new mesh: copy existing patch points
    labelList copiedPatchPoints;
    if (!addToMesh_)
    {
        copiedPatchPoints.setSize(firstLayerDisp.size());
        forAll(firstLayerDisp, patchPointI)
        {
            if (addedPoints_[patchPointI].size())
            {
                label meshPointI = meshPoints[patchPointI];
                label zoneI = mesh_.pointZones().whichZone(meshPointI);
                copiedPatchPoints[patchPointI] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        mesh_.points()[meshPointI],         // point
                        -1,         // master point
                        zoneI,      // zone for point
                        true        // supports a cell
                    )
                );
            }
        }
    }


    // Create points for additional layers
    forAll(firstLayerDisp, patchPointI)
    {
        if (addedPoints_[patchPointI].size())
        {
            label meshPointI = meshPoints[patchPointI];

            label zoneI = mesh_.pointZones().whichZone(meshPointI);

            point pt = mesh_.points()[meshPointI];

            vector disp = firstLayerDisp[patchPointI];

            forAll(addedPoints_[patchPointI], i)
            {
                pt += disp;

                label addedVertI = meshMod.setAction
                (
                    polyAddPoint
                    (
                        pt,         // point
                        (addToMesh_ ? meshPointI : -1), // master point
                        zoneI,      // zone for point
                        true        // supports a cell
                    )
                );

                addedPoints_[patchPointI][i] = addedVertI;

                disp *= expansionRatio[patchPointI];
            }
        }
    }


    //
    // Add cells to all boundaryFaces
    //

    labelListList addedCells(pp.size());

    forAll(pp, patchFaceI)
    {
        if (nFaceLayers[patchFaceI] > 0)
        {
            addedCells[patchFaceI].setSize(nFaceLayers[patchFaceI]);

            label meshFaceI = pp.addressing()[patchFaceI];

            label ownZoneI = mesh_.cellZones().whichZone
            (
                mesh_.faceOwner()[meshFaceI]
            );

            for (label i = 0; i < nFaceLayers[patchFaceI]; i++)
            {
                // Note: add from cell (owner of patch face) or from face?
                // for now add from cell so we can map easily.
                addedCells[patchFaceI][i] = meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,             // master point
                        -1,             // master edge
                        -1,             // master face
                        (addToMesh_ ? mesh_.faceOwner()[meshFaceI] : -1),
                                        //master
                        ownZoneI        // zone for cell
                    )
                );
            }
        }
    }



    // Create faces on top of the original patch faces.
    // These faces are created from original patch faces outwards so in order
    // of increasing cell number. So orientation should be same as original
    // patch face for them to have owner<neighbour.

    layerFaces_.setSize(pp.size());

    forAll(pp.localFaces(), patchFaceI)
    {
        label meshFaceI = pp.addressing()[patchFaceI];

        if (addedCells[patchFaceI].size())
        {
            layerFaces_[patchFaceI].setSize(addedCells[patchFaceI].size() + 1);

            // Get duplicated vertices on the patch face.
            const face& f = pp.localFaces()[patchFaceI];

            face newFace(f.size());

            forAll(addedCells[patchFaceI], i)
            {
                forAll(f, fp)
                {
                    if (addedPoints_[f[fp]].empty())
                    {
                        // Keep original point
                        newFace[fp] =
                        (
                            addToMesh_
                          ? meshPoints[f[fp]]
                          : copiedPatchPoints[f[fp]]
                        );
                    }
                    else
                    {
                        // Get new outside point
                        label offset =
                            addedPoints_[f[fp]].size()
                          - addedCells[patchFaceI].size();
                        newFace[fp] = addedPoints_[f[fp]][i+offset];
                    }
                }


                // Get new neighbour
                label nei;
                label patchI;
                label zoneI = -1;
                bool flip = false;


                if (i == addedCells[patchFaceI].size()-1)
                {
                    // Top layer so is patch face.
                    nei = -1;
                    patchI = patchID[patchFaceI];
                    zoneI = mesh_.faceZones().whichZone(meshFaceI);
                    if (zoneI != -1)
                    {
                        const faceZone& fz = mesh_.faceZones()[zoneI];
                        flip = fz.flipMap()[fz.whichFace(meshFaceI)];
                    }
                }
                else
                {
                    // Internal face between layer i and i+1
                    nei = addedCells[patchFaceI][i+1];
                    patchI = -1;
                }


                layerFaces_[patchFaceI][i+1] = meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,                    // face
                        addedCells[patchFaceI][i],  // owner
                        nei,                        // neighbour
                        -1,                         // master point
                        -1,                         // master edge
                        (addToMesh_ ? meshFaceI : -1), // master face
                        false,                      // flux flip
                        patchI,                     // patch for face
                        zoneI,                      // zone for face
                        flip                        // face zone flip
                    )
                );
            }
        }
    }

    //
    // Modify old patch faces to be on the inside
    //

    if (addToMesh_)
    {
        forAll(pp, patchFaceI)
        {
            if (addedCells[patchFaceI].size())
            {
                label meshFaceI = pp.addressing()[patchFaceI];

                layerFaces_[patchFaceI][0] = meshFaceI;

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        pp[patchFaceI],                 // modified face
                        meshFaceI,                      // label of face
                        mesh_.faceOwner()[meshFaceI],   // owner
                        addedCells[patchFaceI][0],      // neighbour
                        false,                          // face flip
                        -1,                             // patch for face
                        true, //false,                  // remove from zone
                        -1, //zoneI,                    // zone for face
                        false                           // face flip in zone
                    )
                );
            }
        }
    }
    else
    {
        // If creating new mesh: reverse original faces and put them
        // in the exposed patch ID.
        forAll(pp, patchFaceI)
        {
            if (nFaceLayers[patchFaceI] > 0)
            {
                label meshFaceI = pp.addressing()[patchFaceI];
                label zoneI = mesh_.faceZones().whichZone(meshFaceI);
                bool zoneFlip = false;
                if (zoneI != -1)
                {
                    const faceZone& fz = mesh_.faceZones()[zoneI];
                    zoneFlip = !fz.flipMap()[fz.whichFace(meshFaceI)];
                }

                // Reverse and renumber old patch face.
                face f(pp.localFaces()[patchFaceI].reverseFace());
                forAll(f, fp)
                {
                    f[fp] = copiedPatchPoints[f[fp]];
                }

                layerFaces_[patchFaceI][0] = meshMod.setAction
                (
                    polyAddFace
                    (
                        f,                          // modified face
                        addedCells[patchFaceI][0],  // owner
                        -1,                         // neighbour
                        -1,                         // masterPoint
                        -1,                         // masterEdge
                        -1,                         // masterFace
                        true,                       // face flip
                        exposedPatchID[patchFaceI], // patch for face
                        zoneI,                      // zone for face
                        zoneFlip                    // face flip in zone
                    )
                );
            }
        }
    }



    //
    // Create 'side' faces, one per edge that is being extended.
    //

    const labelListList& faceEdges = pp.faceEdges();
    const faceList& localFaces = pp.localFaces();
    const edgeList& edges = pp.edges();

    // Get number of layers per edge. This is 0 if edge is not extruded;
    // max of connected faces otherwise.
    labelList edgeLayers(pp.nEdges());

    {
        // Use list over mesh.nEdges() since syncTools does not yet support
        // partial list synchronisation.
        labelList meshEdgeLayers(mesh_.nEdges(), -1);

        forAll(meshEdges, edgeI)
        {
            const edge& e = edges[edgeI];

            label meshEdgeI = meshEdges[edgeI];

            if ((nPointLayers[e[0]] == 0) && (nPointLayers[e[1]] == 0))
            {
                meshEdgeLayers[meshEdgeI] = 0;
            }
            else
            {
                const labelList& eFaces = pp.edgeFaces()[edgeI];

                forAll(eFaces, i)
                {
                    meshEdgeLayers[meshEdgeI] = max
                    (
                        nFaceLayers[eFaces[i]],
                        meshEdgeLayers[meshEdgeI]
                    );
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh_,
            meshEdgeLayers,
            maxEqOp<label>(),
            label(0)            // initial value
        );

        forAll(meshEdges, edgeI)
        {
            edgeLayers[edgeI] = meshEdgeLayers[meshEdges[edgeI]];
        }
    }


    // Mark off which edges have been extruded
    boolList doneEdge(pp.nEdges(), false);


    // Create faces. Per face walk connected edges and find string of edges
    // between the same two faces and extrude string into a single face.
    forAll(pp, patchFaceI)
    {
        const labelList& fEdges = faceEdges[patchFaceI];

        forAll(fEdges, fp)
        {
            // Get string of edges that needs to be extruded as a single face.
            // Returned as indices in fEdges.
            labelPair indexPair
            (
                getEdgeString
                (
                    pp,
                    globalEdgeFaces,
                    doneEdge,
                    patchFaceI,
                    globalFaces.toGlobal(pp.addressing()[patchFaceI])
                )
            );

            //Pout<< "Found unextruded edges in edges:" << fEdges
            //    << " start:" << indexPair[0]
            //    << " end:" << indexPair[1]
            //    << endl;

            const label startFp = indexPair[0];
            const label endFp = indexPair[1];

            if (startFp != -1)
            {
                // Extrude edges from indexPair[0] up to indexPair[1]
                // (note indexPair = indices of edges. There is one more vertex
                //  than edges)
                const face& f = localFaces[patchFaceI];

                labelList stringedVerts;
                if (endFp >= startFp)
                {
                    stringedVerts.setSize(endFp-startFp+2);
                }
                else
                {
                    stringedVerts.setSize(endFp+f.size()-startFp+2);
                }

                label fp = startFp;

                for (label i = 0; i < stringedVerts.size()-1; i++)
                {
                    stringedVerts[i] = f[fp];
                    doneEdge[fEdges[fp]] = true;
                    fp = f.fcIndex(fp);
                }
                stringedVerts.last() = f[fp];


                // Now stringedVerts contains the vertices in order of face f.
                // This is consistent with the order if f becomes the owner cell
                // and nbrFaceI the neighbour cell. Note that the cells get
                // added in order of pp so we can just use face ordering and
                // because we loop in incrementing order as well we will
                // always have nbrFaceI > patchFaceI.

                label startEdgeI = fEdges[startFp];

                label meshEdgeI = meshEdges[startEdgeI];

                label numEdgeSideFaces = edgeLayers[startEdgeI];

                for (label i = 0; i < numEdgeSideFaces; i++)
                {
                    label vEnd = stringedVerts.last();
                    label vStart = stringedVerts[0];

                    // calculate number of points making up a face
                    label newFp = 2*stringedVerts.size();

                    if (i == 0)
                    {
                        // layer 0 gets all the truncation of neighbouring
                        // faces with more layers.
                        if (addedPoints_[vEnd].size())
                        {
                            newFp +=
                                addedPoints_[vEnd].size() - numEdgeSideFaces;
                        }
                        if (addedPoints_[vStart].size())
                        {
                            newFp +=
                                addedPoints_[vStart].size() - numEdgeSideFaces;
                        }
                    }

                    face newFace(newFp);

                    newFp = 0;

                    // For layer 0 get pp points, for all other layers get
                    // points of layer-1.
                    if (i == 0)
                    {
                        forAll(stringedVerts, stringedI)
                        {
                            label v = stringedVerts[stringedI];
                            addVertex
                            (
                                (
                                    addToMesh_
                                  ? meshPoints[v]
                                  : copiedPatchPoints[v]
                                ),
                                newFace,
                                newFp
                            );
                        }
                    }
                    else
                    {
                        forAll(stringedVerts, stringedI)
                        {
                            label v = stringedVerts[stringedI];
                            if (addedPoints_[v].size())
                            {
                                label offset =
                                    addedPoints_[v].size() - numEdgeSideFaces;
                                addVertex
                                (
                                    addedPoints_[v][i+offset-1],
                                    newFace,
                                    newFp
                                );
                            }
                            else
                            {
                                addVertex
                                (
                                    (
                                        addToMesh_
                                      ? meshPoints[v]
                                      : copiedPatchPoints[v]
                                    ),
                                    newFace,
                                    newFp
                                );
                            }
                        }
                    }

                    // add points between stringed vertices (end)
                    if (numEdgeSideFaces < addedPoints_[vEnd].size())
                    {
                        if (i == 0 && addedPoints_[vEnd].size())
                        {
                            label offset =
                                addedPoints_[vEnd].size() - numEdgeSideFaces;
                            for (label ioff = 0; ioff < offset; ioff++)
                            {
                                addVertex
                                (
                                    addedPoints_[vEnd][ioff],
                                    newFace,
                                    newFp
                                );
                            }
                        }
                    }

                    forAllReverse(stringedVerts, stringedI)
                    {
                        label v = stringedVerts[stringedI];
                        if (addedPoints_[v].size())
                        {
                            label offset =
                                addedPoints_[v].size() - numEdgeSideFaces;
                            addVertex
                            (
                                addedPoints_[v][i+offset],
                                newFace,
                                newFp
                            );
                        }
                        else
                        {
                            addVertex
                            (
                                (
                                    addToMesh_
                                  ? meshPoints[v]
                                  : copiedPatchPoints[v]
                                ),
                                newFace,
                                newFp
                            );
                        }
                    }


                    // add points between stringed vertices (start)
                    if (numEdgeSideFaces < addedPoints_[vStart].size())
                    {
                        if (i == 0 && addedPoints_[vStart].size())
                        {
                            label offset =
                                addedPoints_[vStart].size() - numEdgeSideFaces;
                            for (label ioff = offset-1; ioff >= 0; ioff--)
                            {
                                addVertex
                                (
                                    addedPoints_[vStart][ioff],
                                    newFace,
                                    newFp
                                );
                            }
                        }
                    }

                    if (newFp >= 3)
                    {
                        // Add face inbetween faces patchFaceI and nbrFaceI
                        // (possibly -1 for external edges)

                        newFace.setSize(newFp);

                        if (debug)
                        {
                            labelHashSet verts(2*newFace.size());
                            forAll(newFace, fp)
                            {
                                if (!verts.insert(newFace[fp]))
                                {
                                    FatalErrorIn
                                    (
                                        "addPatchCellLayer::setRefinement(..)"
                                    )   << "Duplicate vertex in face"
                                        << " to be added." << nl
                                        << "newFace:" << newFace << nl
                                        << "points:"
                                        <<  UIndirectList<point>
                                            (
                                                meshMod.points(),
                                                newFace
                                            ) << nl
                                        << "Layer:" << i
                                        << " out of:" << numEdgeSideFaces << nl
                                        << "ExtrudeEdge:" << meshEdgeI
                                        << " at:"
                                        <<  mesh_.edges()[meshEdgeI].line
                                            (
                                                mesh_.points()
                                            ) << nl
                                        << "string:" << stringedVerts
                                        << "stringpoints:"
                                        << UIndirectList<point>
                                            (
                                                pp.localPoints(),
                                                stringedVerts
                                            ) << nl
                                        << "stringNLayers:"
                                        <<  UIndirectList<label>
                                            (
                                                nPointLayers,
                                                stringedVerts
                                            ) << nl
                                        << abort(FatalError);
                                }
                            }
                        }

                        label nbrFaceI = nbrFace
                        (
                            pp.edgeFaces(),
                            startEdgeI,
                            patchFaceI
                        );

                        const labelList& meshFaces = mesh_.edgeFaces
                        (
                            meshEdgeI,
                            ef
                        );

                        // Because we walk in order of patch face and in order
                        // of face edges so face orientation will be opposite
                        // that of the patch edge
                        bool zoneFlip = false;
                        if (edgeZoneID[startEdgeI] != -1)
                        {
                            zoneFlip = !edgeFlip[startEdgeI];
                        }

                        addSideFace
                        (
                            pp,
                            addedCells,

                            newFace,                // vertices of new face
                            edgePatchID[startEdgeI],// -1 or patch for face
                            edgeZoneID[startEdgeI],
                            zoneFlip,
                            inflateFaceID[startEdgeI],

                            patchFaceI,
                            nbrFaceI,
                            meshEdgeI,          // (mesh) edge to inflate
                            i,                  // layer
                            numEdgeSideFaces,   // num layers
                            meshFaces,          // edgeFaces
                            meshMod
                        );
                    }
                }
            }
        }
    }
}


void Foam::addPatchCellLayer::updateMesh
(
    const mapPolyMesh& morphMap,
    const labelList& faceMap,   // new to old patch faces
    const labelList& pointMap   // new to old patch points
)
{
    {
        labelListList newAddedPoints(pointMap.size());

        forAll(newAddedPoints, newPointI)
        {
            label oldPointI = pointMap[newPointI];

            const labelList& added = addedPoints_[oldPointI];

            labelList& newAdded = newAddedPoints[newPointI];
            newAdded.setSize(added.size());
            label newI = 0;

            forAll(added, i)
            {
                label newPointI = morphMap.reversePointMap()[added[i]];

                if (newPointI >= 0)
                {
                    newAdded[newI++] = newPointI;
                }
            }
            newAdded.setSize(newI);
        }
        addedPoints_.transfer(newAddedPoints);
    }

    {
        labelListList newLayerFaces(faceMap.size());

        forAll(newLayerFaces, newFaceI)
        {
            label oldFaceI = faceMap[newFaceI];

            const labelList& added = layerFaces_[oldFaceI];

            labelList& newAdded = newLayerFaces[newFaceI];
            newAdded.setSize(added.size());
            label newI = 0;

            forAll(added, i)
            {
                label newFaceI = morphMap.reverseFaceMap()[added[i]];

                if (newFaceI >= 0)
                {
                    newAdded[newI++] = newFaceI;
                }
            }
            newAdded.setSize(newI);
        }
        layerFaces_.transfer(newLayerFaces);
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2017,2020-2021 OpenCFD Ltd.
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
    const label edgei,
    const label facei
)
{
    const labelList& eFaces = edgeFaces[edgei];

    if (eFaces.size() == 2)
    {
        return (eFaces[0] != facei ? eFaces[0] : eFaces[1]);
    }
    else
    {
        return -1;
    }
}


void Foam::addPatchCellLayer::addVertex
(
    const label pointi,
    face& f,
    label& fp
)
{
    if (fp == 0)
    {
        f[fp++] = pointi;
    }
    else
    {
        if (f[fp-1] != pointi && f[0] != pointi)
        {
            f[fp++] = pointi;
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
    const label thisGlobalFacei,
    const label nbrGlobalFacei,
    const label edgei
) const
{
    const edge& e = pp.edges()[edgei];

    return
        !doneEdge[edgei]                            // not yet handled
     && (
            addedPoints_[e[0]].size()               // is extruded
         || addedPoints_[e[1]].size()
        )
     && (
            nbrFace(globalEdgeFaces, edgei, thisGlobalFacei)
         == nbrGlobalFacei  // is to same neighbour
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
    const label patchFacei,
    const label globalFacei
) const
{
    const labelList& fEdges = pp.faceEdges()[patchFacei];

    label startFp = -1;
    label endFp = -1;

    // Get edge that hasn't been done yet but needs extrusion
    forAll(fEdges, fp)
    {
        label edgei = fEdges[fp];
        const edge& e = pp.edges()[edgei];

        if
        (
            !doneEdge[edgei]
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
        label nbrGlobalFacei = nbrFace
        (
            globalEdgeFaces,
            fEdges[startFp],
            globalFacei
        );

        if (nbrGlobalFacei == -1)
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

            const label initFp = startFp;
            while (true)
            {
                label prevFp = fEdges.rcIndex(startFp);

                if (prevFp == initFp)
                {
                    const edge& e = pp.edges()[fEdges[initFp]];
                    const face& localF = pp.localFaces()[patchFacei];

                    FatalErrorInFunction
                        << "On face:" << patchFacei
                        << " fc:" << pp.faceCentres()[patchFacei]
                        << " vertices:" << localF
                        << " points:"
                        << UIndirectList<point>(pp.points(), pp[patchFacei])
                        << " edges:" << fEdges
                        << " All edges of face seem to have same neighbour "
                        << nbrGlobalFacei
                        << " starting walking from edge " << e
                        << exit(FatalError);
                }

                if
                (
                    !sameEdgeNeighbour
                    (
                        pp,
                        globalEdgeFaces,
                        doneEdge,
                        globalFacei,
                        nbrGlobalFacei,
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
                        globalFacei,
                        nbrGlobalFacei,
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
    const label zonei,
    const bool edgeFlip,
    const label inflateFacei,

    const label ownFacei,               // pp face that provides owner
    const label nbrFacei,
    const label meshEdgei,              // corresponding mesh edge
    const label layeri,                 // layer
    const label numEdgeFaces,           // number of layers for edge
    const labelList& meshFaces,         // precalculated edgeFaces
    polyTopoChange& meshMod
) const
{
    // Adds a side face i.e. extrudes a patch edge.

    label addedFacei = -1;


    // Is patch edge external edge of indirectPrimitivePatch?
    if (nbrFacei == -1)
    {
        // External edge so external face.

        // Determine if different number of layer on owner and neighbour side
        // (relevant only for coupled faces). See section for internal edge
        // below.

        label layerOwn;

        if (addedCells[ownFacei].size() < numEdgeFaces)
        {
            label offset = numEdgeFaces - addedCells[ownFacei].size();
            if (layeri <= offset)
            {
                layerOwn = 0;
            }
            else
            {
                layerOwn = layeri - offset;
            }
        }
        else
        {
            layerOwn = layeri;
        }


        //Pout<< "Added boundary face:" << newFace
        //    << " own:" << addedCells[ownFacei][layerOwn]
        //    << " patch:" << newPatchID
        //    << endl;

        addedFacei = meshMod.setAction
        (
            polyAddFace
            (
                newFace,                        // face
                addedCells[ownFacei][layerOwn], // owner
                -1,                             // neighbour
                -1,                             // master point
                -1,                             // master edge
                inflateFacei,                   // master face
                false,                          // flux flip
                newPatchID,                     // patch for face
                zonei,                          // zone for face
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

        if (addedCells[ownFacei].size() > addedCells[nbrFacei].size())
        {
            label offset =
                addedCells[ownFacei].size() - addedCells[nbrFacei].size();

            layerOwn = layeri;

            if (layeri <= offset)
            {
                layerNbr = 0;
            }
            else
            {
                layerNbr = layeri - offset;
            }
        }
        else if (addedCells[nbrFacei].size() > addedCells[ownFacei].size())
        {
            label offset =
                addedCells[nbrFacei].size() - addedCells[ownFacei].size();

            layerNbr = layeri;

            if (layeri <= offset)
            {
                layerOwn = 0;
            }
            else
            {
                layerOwn = layeri - offset;
            }
        }
        else
        {
            // Same number of layers on both sides.
            layerNbr = layeri;
            layerOwn = layeri;
        }


        // Check mesh internal faces using edge to initialise
        label inflateEdgei = -1;
        if (addToMesh_)
        {
            forAll(meshFaces, i)
            {
                if (mesh_.isInternalFace(meshFaces[i]))
                {
                    // meshEdge uses internal faces so ok to inflate from it
                    inflateEdgei = meshEdgei;
                    break;
                }
            }
        }


        addedFacei = meshMod.setAction
        (
            polyAddFace
            (
                newFace,                        // face
                addedCells[ownFacei][layerOwn], // owner
                addedCells[nbrFacei][layerNbr], // neighbour
                -1,                             // master point
                inflateEdgei,                   // master edge
                -1,                             // master face
                false,                          // flux flip
                -1,                             // patch for face
                zonei,                          // zone for face
                edgeFlip                        // face zone flip
            )
        );

        //Pout<< "Added internal face:" << newFace
        //    << " own:" << addedCells[ownFacei][layerOwn]
        //    << " nei:" << addedCells[nbrFacei][layerNbr]
        //    << endl;
    }

    return addedFacei;
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
        label patchi = mesh.globalData().processorPatches()[i];

        if
        (
            refCast<const processorPolyPatch>(patches[patchi]).neighbProcNo()
         == nbrProcID
        )
        {
            return patchi;
        }
    }
    return -1;
}


void Foam::addPatchCellLayer::setFaceProps
(
    const polyMesh& mesh,
    const label facei,

    label& patchi,
    label& zonei,
    bool& zoneFlip
)
{
    patchi = mesh.boundaryMesh().whichPatch(facei);
    zonei = mesh.faceZones().whichZone(facei);
    if (zonei != -1)
    {
        label index = mesh.faceZones()[zonei].whichFace(facei);
        zoneFlip = mesh.faceZones()[zonei].flipMap()[index];
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
            //FatalErrorInFunction
            WarningInFunction
                << "Problem: cannot find patch edge " << ppEdgeI
                << " with mesh vertices " << patchEdge
                << " at " << patchEdge.line(mesh.points())
                << " in face " << faceI << " with mesh vertices "
                << f
                << " at " << pointField(mesh.points(), f)
                << endl
                << "Continuing with potentially incorrect faceZone orientation"
                //<< exit(FatalError);
                << endl;
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
    const labelUIndList& excludeFaces,
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
            !excludeFaces.found(faceI)
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

    forAll(layerFaces, patchFacei)
    {
        const labelList& faceLabels = layerFaces[patchFacei];

        if (faceLabels.size())
        {
            labelList& added = layerCells[patchFacei];
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
        ListOps::uniqueEqOp<label>(),
        labelList()             // null value
    );

    // Extract pp part
    return labelListList(UIndirectList<labelList>(globalEdgeFaces, meshEdges));
}


void Foam::addPatchCellLayer::markPatchEdges
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    const labelListList& edgeGlobalFaces,
    const labelList& meshEdges,

    bitSet& isPatchEdge,
    bitSet& isPatchBoundaryEdge
)
{
    // Mark (mesh) edges:
    // - anywhere on extrusion
    // - where the extrusion ends

    isPatchEdge.setSize(mesh.nEdges());
    isPatchEdge = false;
    isPatchEdge.set(meshEdges);
    // Synchronise across coupled edges
    syncTools::syncEdgeList
    (
        mesh,
        isPatchEdge,
        orEqOp<unsigned int>(),
        false                   // initial value
    );

    isPatchBoundaryEdge.setSize(mesh.nEdges());
    isPatchBoundaryEdge = false;
    forAll(edgeGlobalFaces, edgei)
    {
        // Test that edge has single global extruded face.
        // Mark on processor that holds the face (since edgeGlobalFaces
        // only gets filled from pp faces so if there is only one this
        // is it)
        if (edgeGlobalFaces[edgei].size() == 1)
        {
            isPatchBoundaryEdge.set(meshEdges[edgei]);
        }
    }
    // Synchronise across coupled edges
    syncTools::syncEdgeList
    (
        mesh,
        isPatchBoundaryEdge,
        orEqOp<unsigned int>(),
        false                   // initial value
    );
}


void Foam::addPatchCellLayer::globalEdgeInfo
(
    const bool zoneFromAnyFace,

    const polyMesh& mesh,
    const globalIndex& globalFaces,
    const labelListList& edgeGlobalFaces,
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,

    labelList& patchEdgeToFace,     // face (in globalFaces index)
    labelList& patchEdgeToPatch,    // patch on face (or -1 for internal faces)
    labelList& patchEdgeToZone,     // zone on face
    bitSet& patchEdgeToFlip         // flip orientation on face
)
{
    // For every edge on the outside of the patch return a potential patch/
    // faceZone to extrude into.

    // Mark (mesh) edges on pp.
    bitSet isExtrudeEdge;
    bitSet isBoundaryEdge;
    markPatchEdges
    (
        mesh,
        pp,
        edgeGlobalFaces,
        meshEdges,

        isExtrudeEdge,
        isBoundaryEdge
    );

    // Build map of pp edges (in mesh point indexing). Note that this
    // is now also on processors that do not have pp (but do have the edge)
    EdgeMap<label> isBoundaryEdgeSet(pp.nEdges());
    for (const label edgei : isBoundaryEdge)
    {
        isBoundaryEdgeSet.insert(mesh.edges()[edgei], edgei);
    }
    EdgeMap<label> isExtrudeEdgeSet(pp.nEdges());
    for (const label edgei : isExtrudeEdge)
    {
        isExtrudeEdgeSet.insert(mesh.edges()[edgei], edgei);
    }


    const faceZoneMesh& fzs = mesh.faceZones();

    // Extract zone info into mesh face indexing for ease of addressing
    labelList faceToZone(mesh.nFaces(), -1);
    bitSet faceToFlip(mesh.nFaces());
    for (const faceZone& fz: fzs)
    {
        const labelList& addressing = fz;
        UIndirectList<label>(faceToZone, addressing) = fz.index();

        const boolList& fm = fz.flipMap();
        forAll(addressing, i)
        {
            faceToFlip[addressing[i]] = fm[i];
        }
    }


    // Storage (over all mesh edges)
    // - face that data originates from (in globalFaces indexing)
    labelList meshEdgeToFace(mesh.nEdges(), -1);
    // - patch (for boundary faces)
    labelList meshEdgeToPatch(mesh.nEdges(), -1);
    // - faceZone
    labelList meshEdgeToZone(mesh.nEdges(), -1);
    // - faceZone orientation
    bitSet meshEdgeToFlip(mesh.nEdges());

    //if (useInternalFaces)
    {
        const bitSet isInternalOrCoupled
        (
            syncTools::getInternalOrCoupledFaces(mesh)
        );

        // Loop over edges of the face to find any faceZone.
        // Edges kept as point pair so we don't construct mesh.faceEdges etc.

        for (const label facei : isInternalOrCoupled)
        {
            const face& f = mesh.faces()[facei];

            label prevPointi = f.last();
            for (const label pointi : f)
            {
                const edge e(prevPointi, pointi);

                // Check if edge is internal to extrusion. Take over faceZone
                // etc from internal face.
                const auto eFnd = isExtrudeEdgeSet.cfind(e);
                if (eFnd)
                {
                    const label edgei = eFnd();

                    if (faceToZone[facei] != -1)
                    {
                        // Found a zoned internal face. Use.
                        meshEdgeToFace[edgei] = globalFaces.toGlobal(facei);
                        meshEdgeToZone[edgei] = faceToZone[facei];
                        const edge& meshE = mesh.edges()[edgei];
                        const int d = edge::compare(e, meshE);
                        if (d == 1)
                        {
                            meshEdgeToFlip[edgei] = faceToFlip[facei];
                        }
                        else if (d == -1)
                        {
                            meshEdgeToFlip[edgei] = !faceToFlip[facei];
                        }
                        else
                        {
                            FatalErrorInFunction << "big problem"
                                << exit(FatalError);
                        }
                    }
                }

                prevPointi = pointi;
            }
        }
    }


    //if (useBoundaryFaces)
    {
        // Loop over all patches and find 'best' one (non-coupled,
        // non-extrusion, non-constraint(?)). Note that logic is slightly
        // different from internal faces loop above - first patch face
        // is being used instead of first zoned face for internal faces

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        bitSet isPpFace(mesh.nFaces());
        isPpFace.set(pp.addressing());
        // Note: no need to sync ppFace since does not include processor patches

        for (const polyPatch& pp : patches)
        {
            if (!pp.coupled())
            {
                // TBD. Check for constraint? This is usually a geometric
                //      constraint and not a topological one so should
                //      be handled in the extrusion vector calculation instead?

                forAll(pp, i)
                {
                    const label facei = pp.start()+i;

                    if (!isPpFace[facei])
                    {
                        const face& f = pp[i];

                        label prevPointi = f.last();
                        for (const label pointi : f)
                        {
                            const edge e(prevPointi, pointi);
                            const auto eFnd =
                            (
                                zoneFromAnyFace
                              ? isExtrudeEdgeSet.cfind(e)
                              : isBoundaryEdgeSet.cfind(e)
                            );
                            if (eFnd)
                            {
                                const label edgei = eFnd();
                                if (meshEdgeToFace[edgei] == -1)
                                {
                                    // Found unassigned face. Use its
                                    // information.
                                    // Note that we use the lowest numbered
                                    // patch face.

                                    meshEdgeToFace[edgei] =
                                        globalFaces.toGlobal(facei);

                                    // Override any patch info
                                    if (meshEdgeToPatch[edgei] == -1)
                                    {
                                        meshEdgeToPatch[edgei] = pp.index();
                                    }

                                    // Override any zone info
                                    if (meshEdgeToZone[edgei] == -1)
                                    {
                                        meshEdgeToZone[edgei] =
                                            faceToZone[facei];
                                        const edge& meshE = mesh.edges()[edgei];
                                        const int d = edge::compare(e, meshE);
                                        if (d == 1)
                                        {
                                            meshEdgeToFlip[edgei] =
                                                faceToFlip[facei];
                                        }
                                        else if (d == -1)
                                        {
                                            meshEdgeToFlip[edgei] =
                                               !faceToFlip[facei];
                                        }
                                        else
                                        {
                                            FatalErrorInFunction
                                                << "big problem"
                                                << exit(FatalError);
                                        }
                                    }
                                }
                            }

                            prevPointi = pointi;
                        }
                    }
                }
            }
        }
    }

    // Synchronise across coupled edges. Max patch/face/faceZone wins
    syncTools::syncEdgeList
    (
        mesh,
        meshEdgeToFace,
        maxEqOp<label>(),
        label(-1)
    );
    syncTools::syncEdgeList
    (
        mesh,
        meshEdgeToPatch,
        maxEqOp<label>(),
        label(-1)
    );
    syncTools::syncEdgeList
    (
        mesh,
        meshEdgeToZone,
        maxEqOp<label>(),
        label(-1)
    );
//    // Note: flipMap not yet done. Needs edge orientation. This is handled
//    // later on.
//    if (Pstream::parRun())
//    {
//        const globalMeshData& gd = mesh.globalData();
//        const indirectPrimitivePatch& cpp = gd.coupledPatch();
//
//        labelList patchEdges;
//        labelList coupledEdges;
//        bitSet sameEdgeOrientation;
//        PatchTools::matchEdges
//        (
//            pp,
//            cpp,
//            patchEdges,
//            coupledEdges,
//            sameEdgeOrientation
//        );
//
//        // Convert data on pp edges to data on coupled patch
//        labelList cppEdgeZoneID(cpp.nEdges(), -1);
//        boolList cppEdgeFlip(cpp.nEdges(), false);
//        forAll(coupledEdges, i)
//        {
//            label cppEdgei = coupledEdges[i];
//            label ppEdgei = patchEdges[i];
//
//            cppEdgeZoneID[cppEdgei] = edgeZoneID[ppEdgei];
//            if (sameEdgeOrientation[i])
//            {
//                cppEdgeFlip[cppEdgei] = edgeFlip[ppEdgei];
//            }
//            else
//            {
//                cppEdgeFlip[cppEdgei] = !edgeFlip[ppEdgei];
//            }
//        }
//
//        // Sync
//        const globalIndexAndTransform& git = gd.globalTransforms();
//        const mapDistribute& edgeMap = gd.globalEdgeSlavesMap();
//
//        globalMeshData::syncData
//        (
//            cppEdgeZoneID,
//            gd.globalEdgeSlaves(),
//            gd.globalEdgeTransformedSlaves(),
//            edgeMap,
//            git,
//            maxEqOp<label>(),
//            dummyTransform()
//        );
//        globalMeshData::syncData
//        (
//            cppEdgeFlip,
//            gd.globalEdgeSlaves(),
//            gd.globalEdgeTransformedSlaves(),
//            edgeMap,
//            git,
//            andEqOp<bool>(),
//            dummyTransform()
//        );
//
//        // Convert data on coupled edges to pp edges
//        forAll(coupledEdges, i)
//        {
//            label cppEdgei = coupledEdges[i];
//            label ppEdgei = patchEdges[i];
//
//            edgeZoneID[ppEdgei] = cppEdgeZoneID[cppEdgei];
//            if (sameEdgeOrientation[i])
//            {
//                edgeFlip[ppEdgei] = cppEdgeFlip[cppEdgei];
//            }
//            else
//            {
//                edgeFlip[ppEdgei] = !cppEdgeFlip[cppEdgei];
//            }
//        }
//    }

    syncTools::syncEdgeList
    (
        mesh,
        meshEdgeToFlip,
        orEqOp<unsigned int>(),
        0
    );

    // Extract pp info
    patchEdgeToFace = UIndirectList<label>(meshEdgeToFace, meshEdges);
    patchEdgeToPatch = UIndirectList<label>(meshEdgeToPatch, meshEdges);
    patchEdgeToZone = UIndirectList<label>(meshEdgeToZone, meshEdges);
    patchEdgeToFlip.setSize(meshEdges.size());
    patchEdgeToFlip = false;
    forAll(meshEdges, i)
    {
        patchEdgeToFlip[i] = meshEdgeToFlip[meshEdges[i]];
    }
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
    nPatches = patches.size();
    edgeZoneID.setSize(pp.nEdges());
    edgeZoneID = -1;
    edgeFlip.setSize(pp.nEdges());
    edgeFlip = false;
    inflateFaceID.setSize(pp.nEdges(), -1);


    // Determine properties for faces extruded from edges
    // - edge inbetween two different processors:
    //      - extrude as patch face on correct processor
    // - edge at end of patch (so edgeFaces.size() == 1):
    //      - use mesh face that is a boundary face
    // - edge internal to patch (so edgeFaces.size() == 2):


    // Pass1:
    // For all edges inbetween two processors: see if matches to existing
    // processor patch or create interprocessor-patch if necessary.
    // Sets edgePatchID[edgeI] but none of the other quantities

    forAll(globalEdgeFaces, edgei)
    {
        const labelList& eGlobalFaces = globalEdgeFaces[edgei];
        if
        (
            eGlobalFaces.size() == 2
         && pp.edgeFaces()[edgei].size() == 1
        )
        {
            // Locally but not globally a boundary edge. Hence a coupled
            // edge. Find the patch to use if on different processors.

            label f0 = eGlobalFaces[0];
            label f1 = eGlobalFaces[1];

            label otherProci = -1;
            if (globalFaces.isLocal(f0) && !globalFaces.isLocal(f1))
            {
                otherProci = globalFaces.whichProcID(f1);
            }
            else if (!globalFaces.isLocal(f0) && globalFaces.isLocal(f1))
            {
                otherProci = globalFaces.whichProcID(f0);
            }


            if (otherProci != -1)
            {
                if (gd[Pstream::myProcNo()].found(otherProci))
                {
                    // There is already a processorPolyPatch to otherProci.
                    // Use it. Note that we can only index procPatchMap
                    // if the processor actually is a neighbour processor
                    // so that is why we first check.
                    edgePatchID[edgei] = gd.procPatchMap()[otherProci];
                }
                else
                {
                    // Cannot find a patch to processor. See if already
                    // marked for addition
                    if (nbrProcToPatch.found(otherProci))
                    {
                        edgePatchID[edgei] = nbrProcToPatch[otherProci];
                    }
                    else
                    {
                        edgePatchID[edgei] = nPatches;
                        nbrProcToPatch.insert(otherProci, nPatches);
                        patchToNbrProc.insert(nPatches, otherProci);
                        nPatches++;
                    }
                }
            }
        }
    }


    // Pass2: determine face properties for all other edges
    // ----------------------------------------------------

    // Info for edges of pp
    labelList edgeToFace;
    labelList edgeToPatch;
    labelList edgeToZone;
    bitSet edgeToFlip;
    globalEdgeInfo
    (
        zoneFromAnyFace,    // internal edge info also from boundary faces

        mesh,
        globalFaces,
        globalEdgeFaces,
        pp,
        meshEdges,

        edgeToFace,         // face (in globalFaces index)
        edgeToPatch,        // patch on face (or -1 for internal faces)
        edgeToZone,         // zone on face
        edgeToFlip          // flip orientation on face
    );

    const labelListList& edgeFaces = pp.edgeFaces();

    DynamicList<label> dynMeshEdgeFaces;

    forAll(edgeFaces, edgei)
    {
        if (edgePatchID[edgei] == -1)
        {
            if (edgeFaces[edgei].size() == 2)
            {
                // Internal edge. Look at any face (internal or boundary) to
                // determine extrusion properties. First one that has zone
                // info wins
                if (globalFaces.isLocal(edgeToFace[edgei]))
                {
                    inflateFaceID[edgei] =
                        globalFaces.toLocal(edgeToFace[edgei]);
                }
                edgeZoneID[edgei] = edgeToZone[edgei];
                edgeFlip[edgei] = edgeToFlip[edgei];
            }
            else
            {
                // Proper, uncoupled patch edge. Boundary to get info from
                // might be on a different processor!

                if (globalFaces.isLocal(edgeToFace[edgei]))
                {
                    inflateFaceID[edgei] =
                        globalFaces.toLocal(edgeToFace[edgei]);
                }
                edgePatchID[edgei] = edgeToPatch[edgei];
                edgeZoneID[edgei] = edgeToZone[edgei];
                edgeFlip[edgei] = edgeToFlip[edgei];
            }
        }
    }



    // Now hopefully every boundary edge has a edge patch. Check
    if (debug)
    {
        forAll(edgeFaces, edgei)
        {
            if (edgeFaces[edgei].size() == 1 && edgePatchID[edgei] == -1)
            {
                const edge& e = pp.edges()[edgei];
                WarningInFunction
                    << "Have no sidePatchID for edge " << edgei << " points "
                    << pp.points()[pp.meshPoints()[e[0]]]
                    << pp.points()[pp.meshPoints()[e[1]]]
                    << endl;
            }
        }
    }



    // Pass3: for any faces set in pass1 see if we can find a processor face
    //        to inherit from (we only have a patch, not a patch face)
    forAll(edgeFaces, edgei)
    {
        if
        (
            edgeFaces[edgei].size() == 1
         && globalEdgeFaces[edgei].size() == 2
         && edgePatchID[edgei] != -1
         && inflateFaceID[edgei] == -1
        )
        {
            // 1. Do we have a local boundary face to inflate from

            label myFaceI = pp.addressing()[edgeFaces[edgei][0]];

            // Pick up any boundary face on this edge and use its properties
            label meshEdgei = meshEdges[edgei];
            const labelList& meshFaces = mesh.edgeFaces
            (
                meshEdgei,
                dynMeshEdgeFaces
            );

            forAll(meshFaces, k)
            {
                label facei = meshFaces[k];

                if (facei != myFaceI && !mesh.isInternalFace(facei))
                {
                    if (patches.whichPatch(facei) == edgePatchID[edgei])
                    {
                        setFaceProps
                        (
                            mesh,
                            pp,
                            edgei,
                            facei,

                            edgePatchID[edgei],
                            edgeZoneID[edgei],
                            edgeFlip[edgei],
                            inflateFaceID[edgei]
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
        bitSet sameEdgeOrientation;
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
            label cppEdgei = coupledEdges[i];
            label ppEdgei = patchEdges[i];

            cppEdgeZoneID[cppEdgei] = edgeZoneID[ppEdgei];
            if (sameEdgeOrientation[i])
            {
                cppEdgeFlip[cppEdgei] = edgeFlip[ppEdgei];
            }
            else
            {
                cppEdgeFlip[cppEdgei] = !edgeFlip[ppEdgei];
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
            label cppEdgei = coupledEdges[i];
            label ppEdgei = patchEdges[i];

            edgeZoneID[ppEdgei] = cppEdgeZoneID[cppEdgei];
            if (sameEdgeOrientation[i])
            {
                edgeFlip[ppEdgei] = cppEdgeFlip[cppEdgei];
            }
            else
            {
                edgeFlip[ppEdgei] = !cppEdgeFlip[cppEdgei];
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
        FatalErrorInFunction
            << "Size of new points is not same as number of points used by"
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
            FatalErrorInFunction
                << "Illegal number of layers " << nPointLayers[i]
                << " at patch point " << i << abort(FatalError);
        }
    }
    forAll(nFaceLayers, i)
    {
        if (nFaceLayers[i] < 0)
        {
            FatalErrorInFunction
                << "Illegal number of layers " << nFaceLayers[i]
                << " at patch face " << i << abort(FatalError);
        }
    }

    forAll(globalEdgeFaces, edgei)
    {
        if (globalEdgeFaces[edgei].size() > 2)
        {
            const edge& e = pp.edges()[edgei];

            if (nPointLayers[e[0]] > 0 || nPointLayers[e[1]] > 0)
            {
                FatalErrorInFunction
                    << "Trying to extrude edge "
                    << e.line(pp.localPoints())
                    << " which is non-manifold (has "
                    << globalEdgeFaces[edgei].size()
                    << " local or coupled faces using it)"
                    << " of which "
                    << pp.edgeFaces()[edgei].size()
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
            labelList n(mesh_.nPoints(), Zero);
            labelUIndList(n, meshPoints) = nPointLayers;
            syncTools::syncPointList(mesh_, n, maxEqOp<label>(), label(0));

            // Non-synced
            forAll(meshPoints, i)
            {
                label meshPointi = meshPoints[i];

                if (n[meshPointi] != nPointLayers[i])
                {
                    FatalErrorInFunction
                        << "At mesh point:" << meshPointi
                        << " coordinate:" << mesh_.points()[meshPointi]
                        << " specified nLayers:" << nPointLayers[i] << endl
                        << "On coupled point a different nLayers:"
                        << n[meshPointi] << " was specified."
                        << abort(FatalError);
                }
            }


            // Check that nPointLayers equals the max layers of connected faces
            // (or 0). Anything else makes no sense.
            labelList nFromFace(mesh_.nPoints(), Zero);
            forAll(nFaceLayers, i)
            {
                const face& f = pp[i];

                forAll(f, fp)
                {
                    label pointi = f[fp];

                    nFromFace[pointi] = max(nFromFace[pointi], nFaceLayers[i]);
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
                label meshPointi = meshPoints[i];

                if
                (
                    nPointLayers[i] > 0
                 && nPointLayers[i] != nFromFace[meshPointi]
                )
                {
                    FatalErrorInFunction
                        << "At mesh point:" << meshPointi
                        << " coordinate:" << mesh_.points()[meshPointi]
                        << " specified nLayers:" << nPointLayers[i] << endl
                        << "but the max nLayers of surrounding faces is:"
                        << nFromFace[meshPointi]
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
                label meshPointi = meshPoints[i];

                if (mag(d[meshPointi] - firstLayerDisp[i]) > SMALL)
                {
                    FatalErrorInFunction
                        << "At mesh point:" << meshPointi
                        << " coordinate:" << mesh_.points()[meshPointi]
                        << " specified displacement:" << firstLayerDisp[i]
                        << endl
                        << "On coupled point a different displacement:"
                        << d[meshPointi] << " was specified."
                        << abort(FatalError);
                }
            }
        }

        // Check that edges of pp (so ones that become boundary faces)
        // connect to only one boundary face. Guarantees uniqueness of
        // patch that they go into so if this is a coupled patch both
        // sides decide the same.
        // ~~~~~~~~~~~~~~~~~~~~~~

        for (label edgei = pp.nInternalEdges(); edgei < pp.nEdges(); edgei++)
        {
            const edge& e = pp.edges()[edgei];

            if (nPointLayers[e[0]] > 0 || nPointLayers[e[1]] > 0)
            {
                // Edge is to become a face

                const labelList& eFaces = pp.edgeFaces()[edgei];

                // First check: pp should be single connected.
                if (eFaces.size() != 1)
                {
                    FatalErrorInFunction
                        << "boundary-edge-to-be-extruded:"
                        << pp.points()[meshPoints[e[0]]]
                        << pp.points()[meshPoints[e[1]]]
                        << " has more than two faces using it:" << eFaces
                        << abort(FatalError);
                }

                label myFacei = pp.addressing()[eFaces[0]];

                label meshEdgei = meshEdges[edgei];

                // Mesh faces using edge
                const labelList& meshFaces = mesh_.edgeFaces(meshEdgei, ef);

                // Check that there is only one patchface using edge.
                const polyBoundaryMesh& patches = mesh_.boundaryMesh();

                label bFacei = -1;

                forAll(meshFaces, i)
                {
                    label facei = meshFaces[i];

                    if (facei != myFacei)
                    {
                        if (!mesh_.isInternalFace(facei))
                        {
                            if (bFacei == -1)
                            {
                                bFacei = facei;
                            }
                            else
                            {
                                FatalErrorInFunction
                                    << "boundary-edge-to-be-extruded:"
                                    << pp.points()[meshPoints[e[0]]]
                                    << pp.points()[meshPoints[e[1]]]
                                    << " has more than two boundary faces"
                                    << " using it:"
                                    << bFacei << " fc:"
                                    << mesh_.faceCentres()[bFacei]
                                    << " patch:" << patches.whichPatch(bFacei)
                                    << " and " << facei << " fc:"
                                    << mesh_.faceCentres()[facei]
                                    << " patch:" << patches.whichPatch(facei)
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

    forAll(pp, patchFacei)
    {
        label meshFacei = pp.addressing()[patchFacei];

        patchID[patchFacei] = patches.whichPatch(meshFacei);
    }


    // From master point (in patch point label) to added points (in mesh point
    // label)
    addedPoints_.setSize(pp.nPoints());

    // Mark points that do not get extruded by setting size of addedPoints_ to 0
    label nTruncated = 0;

    forAll(nPointLayers, patchPointi)
    {
        if (nPointLayers[patchPointi] > 0)
        {
            addedPoints_[patchPointi].setSize(nPointLayers[patchPointi]);
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
        forAll(firstLayerDisp, patchPointi)
        {
            label meshPointi = meshPoints[patchPointi];
            label zoneI = mesh_.pointZones().whichZone(meshPointi);
            copiedPatchPoints[patchPointi] = meshMod.setAction
            (
                polyAddPoint
                (
                    mesh_.points()[meshPointi],         // point
                    -1,         // master point
                    zoneI,      // zone for point
                    true        // supports a cell
                )
            );
        }
    }


    // Create points for additional layers
    forAll(firstLayerDisp, patchPointi)
    {
        if (addedPoints_[patchPointi].size())
        {
            label meshPointi = meshPoints[patchPointi];

            label zoneI = mesh_.pointZones().whichZone(meshPointi);

            point pt = mesh_.points()[meshPointi];

            vector disp = firstLayerDisp[patchPointi];

            forAll(addedPoints_[patchPointi], i)
            {
                pt += disp;

                label addedVertI = meshMod.setAction
                (
                    polyAddPoint
                    (
                        pt,         // point
                        (addToMesh_ ? meshPointi : -1), // master point
                        zoneI,      // zone for point
                        true        // supports a cell
                    )
                );

                addedPoints_[patchPointi][i] = addedVertI;

                disp *= expansionRatio[patchPointi];
            }
        }
    }


    //
    // Add cells to all boundaryFaces
    //

    labelListList addedCells(pp.size());

    forAll(pp, patchFacei)
    {
        if (nFaceLayers[patchFacei] > 0)
        {
            addedCells[patchFacei].setSize(nFaceLayers[patchFacei]);

            label meshFacei = pp.addressing()[patchFacei];

            label ownZoneI = mesh_.cellZones().whichZone
            (
                mesh_.faceOwner()[meshFacei]
            );

            for (label i = 0; i < nFaceLayers[patchFacei]; i++)
            {
                // Note: add from cell (owner of patch face) or from face?
                // for now add from cell so we can map easily.
                addedCells[patchFacei][i] = meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,             // master point
                        -1,             // master edge
                        -1,             // master face
                        (addToMesh_ ? mesh_.faceOwner()[meshFacei] : -1),
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

    forAll(pp.localFaces(), patchFacei)
    {
        label meshFacei = pp.addressing()[patchFacei];

        if (addedCells[patchFacei].size())
        {
            layerFaces_[patchFacei].setSize(addedCells[patchFacei].size() + 1);

            // Get duplicated vertices on the patch face.
            const face& f = pp.localFaces()[patchFacei];

            face newFace(f.size());

            forAll(addedCells[patchFacei], i)
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
                          - addedCells[patchFacei].size();
                        newFace[fp] = addedPoints_[f[fp]][i+offset];
                    }
                }


                // Get new neighbour
                label nei;
                label patchi;
                label zoneI = -1;
                bool flip = false;


                if (i == addedCells[patchFacei].size()-1)
                {
                    // Top layer so is patch face.
                    nei = -1;
                    patchi = patchID[patchFacei];
                    zoneI = mesh_.faceZones().whichZone(meshFacei);
                    if (zoneI != -1)
                    {
                        const faceZone& fz = mesh_.faceZones()[zoneI];
                        flip = fz.flipMap()[fz.whichFace(meshFacei)];
                    }
                }
                else
                {
                    // Internal face between layer i and i+1
                    nei = addedCells[patchFacei][i+1];
                    patchi = -1;
                }


                layerFaces_[patchFacei][i+1] = meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,                    // face
                        addedCells[patchFacei][i],  // owner
                        nei,                        // neighbour
                        -1,                         // master point
                        -1,                         // master edge
                        (addToMesh_ ? meshFacei : -1), // master face
                        false,                      // flux flip
                        patchi,                     // patch for face
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
        forAll(pp, patchFacei)
        {
            if (addedCells[patchFacei].size())
            {
                label meshFacei = pp.addressing()[patchFacei];

                layerFaces_[patchFacei][0] = meshFacei;

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        pp[patchFacei],                 // modified face
                        meshFacei,                      // label of face
                        mesh_.faceOwner()[meshFacei],   // owner
                        addedCells[patchFacei][0],      // neighbour
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
        forAll(pp, patchFacei)
        {
            if (nFaceLayers[patchFacei] > 0)
            {
                label meshFacei = pp.addressing()[patchFacei];
                label zoneI = mesh_.faceZones().whichZone(meshFacei);
                bool zoneFlip = false;
                if (zoneI != -1)
                {
                    const faceZone& fz = mesh_.faceZones()[zoneI];
                    zoneFlip = !fz.flipMap()[fz.whichFace(meshFacei)];
                }

                // Reverse and renumber old patch face.
                face f(pp.localFaces()[patchFacei].reverseFace());
                forAll(f, fp)
                {
                    f[fp] = copiedPatchPoints[f[fp]];
                }

                layerFaces_[patchFacei][0] = meshMod.setAction
                (
                    polyAddFace
                    (
                        f,                          // modified face
                        addedCells[patchFacei][0],  // owner
                        -1,                         // neighbour
                        -1,                         // masterPoint
                        -1,                         // masterEdge
                        -1,                         // masterFace
                        true,                       // face flip
                        exposedPatchID[patchFacei], // patch for face
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

        forAll(meshEdges, edgei)
        {
            const edge& e = edges[edgei];

            label meshEdgei = meshEdges[edgei];

            if ((nPointLayers[e[0]] == 0) && (nPointLayers[e[1]] == 0))
            {
                meshEdgeLayers[meshEdgei] = 0;
            }
            else
            {
                const labelList& eFaces = pp.edgeFaces()[edgei];

                forAll(eFaces, i)
                {
                    meshEdgeLayers[meshEdgei] = max
                    (
                        nFaceLayers[eFaces[i]],
                        meshEdgeLayers[meshEdgei]
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

        forAll(meshEdges, edgei)
        {
            edgeLayers[edgei] = meshEdgeLayers[meshEdges[edgei]];
        }
    }


    // Mark off which edges have been extruded
    boolList doneEdge(pp.nEdges(), false);


    // Create faces. Per face walk connected edges and find string of edges
    // between the same two faces and extrude string into a single face.
    forAll(pp, patchFacei)
    {
        const labelList& fEdges = faceEdges[patchFacei];

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
                    patchFacei,
                    globalFaces.toGlobal(pp.addressing()[patchFacei])
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
                const face& f = localFaces[patchFacei];

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
                // and nbrFacei the neighbour cell. Note that the cells get
                // added in order of pp so we can just use face ordering and
                // because we loop in incrementing order as well we will
                // always have nbrFacei > patchFacei.

                label startEdgei = fEdges[startFp];

                label meshEdgei = meshEdges[startEdgei];

                label numEdgeSideFaces = edgeLayers[startEdgei];

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
                        // Add face inbetween faces patchFacei and nbrFacei
                        // (possibly -1 for external edges)

                        newFace.setSize(newFp);

                        if (debug)
                        {
                            labelHashSet verts(2*newFace.size());
                            forAll(newFace, fp)
                            {
                                if (!verts.insert(newFace[fp]))
                                {
                                    FatalErrorInFunction
                                        << "Duplicate vertex in face"
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
                                        << "ExtrudeEdge:" << meshEdgei
                                        << " at:"
                                        <<  mesh_.edges()[meshEdgei].line
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
                                        <<  labelUIndList
                                            (
                                                nPointLayers,
                                                stringedVerts
                                            ) << nl
                                        << abort(FatalError);
                                }
                            }
                        }

                        label nbrFacei = nbrFace
                        (
                            pp.edgeFaces(),
                            startEdgei,
                            patchFacei
                        );

                        const labelList& meshFaces = mesh_.edgeFaces
                        (
                            meshEdgei,
                            ef
                        );

                        // Because we walk in order of patch face and in order
                        // of face edges so face orientation will be opposite
                        // that of the patch edge
                        bool zoneFlip = false;
                        if (edgeZoneID[startEdgei] != -1)
                        {
                            zoneFlip = !edgeFlip[startEdgei];
                        }

                        addSideFace
                        (
                            pp,
                            addedCells,

                            newFace,                // vertices of new face
                            edgePatchID[startEdgei],// -1 or patch for face
                            edgeZoneID[startEdgei],
                            zoneFlip,
                            inflateFaceID[startEdgei],

                            patchFacei,
                            nbrFacei,
                            meshEdgei,          // (mesh) edge to inflate
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

        forAll(newAddedPoints, newPointi)
        {
            label oldPointi = pointMap[newPointi];

            const labelList& added = addedPoints_[oldPointi];

            labelList& newAdded = newAddedPoints[newPointi];
            newAdded.setSize(added.size());
            label newI = 0;

            forAll(added, i)
            {
                label newPointi = morphMap.reversePointMap()[added[i]];

                if (newPointi >= 0)
                {
                    newAdded[newI++] = newPointi;
                }
            }
            newAdded.setSize(newI);
        }
        addedPoints_.transfer(newAddedPoints);
    }

    {
        labelListList newLayerFaces(faceMap.size());

        forAll(newLayerFaces, newFacei)
        {
            label oldFacei = faceMap[newFacei];

            const labelList& added = layerFaces_[oldFacei];

            labelList& newAdded = newLayerFaces[newFacei];
            newAdded.setSize(added.size());
            label newI = 0;

            forAll(added, i)
            {
                label newFacei = morphMap.reverseFaceMap()[added[i]];

                if (newFacei >= 0)
                {
                    newAdded[newI++] = newFacei;
                }
            }
            newAdded.setSize(newI);
        }
        layerFaces_.transfer(newLayerFaces);
    }
}


// ************************************************************************* //

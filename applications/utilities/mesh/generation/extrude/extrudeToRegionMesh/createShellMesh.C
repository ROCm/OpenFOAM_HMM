/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "createShellMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "polyAddPoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyAddCell.H"
#include "patchPointEdgeCirculator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::createShellMesh, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::createShellMesh::calcPointRegions
(
    const primitiveFacePatch& patch,
    const PackedBoolList& nonManifoldEdge,
    faceList& pointRegions,
    labelList& regionPoints
)
{
    pointRegions.setSize(patch.size());
    forAll(pointRegions, faceI)
    {
        const face& f = patch.localFaces()[faceI];
        pointRegions[faceI].setSize(f.size(), -1);
    }

    label nRegions = 0;

    forAll(pointRegions, faceI)
    {
        const face& f = patch.localFaces()[faceI];

        forAll(f, fp)
        {
            if (pointRegions[faceI][fp] == -1)
            {
                // Found unassigned point. Distribute current region.
                label pointI = f[fp];
                label edgeI = patch.faceEdges()[faceI][fp];

                patchPointEdgeCirculator circ
                (
                    patch,
                    nonManifoldEdge,
                    edgeI,
                    findIndex(patch.edgeFaces()[edgeI], faceI),
                    pointI
                );

                for
                (
                    patchPointEdgeCirculator iter = circ.begin();
                    iter != circ.end();
                    ++iter
                )
                {
                    label face2 = iter.faceID();

                    if (face2 != -1)
                    {
                        const face& f2 = patch.localFaces()[face2];
                        label fp2 = findIndex(f2, pointI);
                        label& region = pointRegions[face2][fp2];
                        if (region != -1)
                        {
                            FatalErrorIn
                            (
                                "createShellMesh::calcPointRegions(..)"
                            )   << "On point " << pointI
                                << " at:" << patch.localPoints()[pointI]
                                << " found region:" << region
                                << abort(FatalError);
                        }
                        region = nRegions;
                    }
                }

                nRegions++;
            }
        }
    }


    // From region back to originating point (many to one, a point might
    // have multiple regions though)
    regionPoints.setSize(nRegions);
    forAll(pointRegions, faceI)
    {
        const face& f = patch.localFaces()[faceI];

        forAll(f, fp)
        {
            regionPoints[pointRegions[faceI][fp]] = f[fp];
        }
    }


    if (debug)
    {
        const labelListList& pointFaces = patch.pointFaces();
        forAll(pointFaces, pointI)
        {
            label region = -1;
            const labelList& pFaces = pointFaces[pointI];
            forAll(pFaces, i)
            {
                label faceI = pFaces[i];
                const face& f = patch.localFaces()[faceI];
                label fp = findIndex(f, pointI);

                if (region == -1)
                {
                    region = pointRegions[faceI][fp];
                }
                else if (region != pointRegions[faceI][fp])
                {
                    Pout<< "Non-manifold point:" << pointI
                        << " at " << patch.localPoints()[pointI]
                        << " region:" << region
                        << " otherRegion:" << pointRegions[faceI][fp]
                        << endl;

                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::createShellMesh::createShellMesh
(
    const primitiveFacePatch& patch,
    const faceList& pointRegions,
    const labelList& regionPoints
)
:
    patch_(patch),
    pointRegions_(pointRegions),
    regionPoints_(regionPoints)
{
    if (pointRegions_.size() != patch_.size())
    {
        FatalErrorIn("createShellMesh::createShellMesh(..)")
            << "nFaces:" << patch_.size()
            << " pointRegions:" << pointRegions.size()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::createShellMesh::setRefinement
(
    const pointField& firstLayerDisp,
    const scalar expansionRatio,
    const label nLayers,
    const labelList& topPatchID,
    const labelList& bottomPatchID,
    const labelListList& extrudeEdgePatches,
    polyTopoChange& meshMod
)
{
    if (firstLayerDisp.size() != regionPoints_.size())
    {
        FatalErrorIn("createShellMesh::setRefinement(..)")
            << "nRegions:" << regionPoints_.size()
            << " firstLayerDisp:" << firstLayerDisp.size()
            << exit(FatalError);
    }

    if
    (
        topPatchID.size() != patch_.size()
     && bottomPatchID.size() != patch_.size()
    )
    {
        FatalErrorIn("createShellMesh::setRefinement(..)")
            << "nFaces:" << patch_.size()
            << " topPatchID:" << topPatchID.size()
            << " bottomPatchID:" << bottomPatchID.size()
            << exit(FatalError);
    }

    if (extrudeEdgePatches.size() != patch_.nEdges())
    {
        FatalErrorIn("createShellMesh::setRefinement(..)")
            << "nEdges:" << patch_.nEdges()
            << " extrudeEdgePatches:" << extrudeEdgePatches.size()
            << exit(FatalError);
    }



    // From cell to patch (trivial)
    DynamicList<label> cellToFaceMap(nLayers*patch_.size());
    // From face to patch+turning index
    DynamicList<label> faceToFaceMap
    (
        (nLayers+1)*(patch_.size()+patch_.nEdges())
    );
    // From face to patch edge index
    DynamicList<label> faceToEdgeMap(nLayers*(patch_.nEdges()+patch_.nEdges()));
    // From point to patch point index
    DynamicList<label> pointToPointMap((nLayers+1)*patch_.nPoints());


    // Introduce new cell for every face
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList addedCells(nLayers*patch_.size());
    forAll(patch_, faceI)
    {
        for (label layerI = 0; layerI < nLayers; layerI++)
        {
            addedCells[nLayers*faceI+layerI] = meshMod.addCell
            (
                -1,                     // masterPointID
                -1,                     // masterEdgeID
                -1,                     // masterFaceID
                cellToFaceMap.size(),   // masterCellID
                -1                      // zoneID
            );
            cellToFaceMap.append(faceI);
        }
    }


    // Introduce original points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Original point numbers in local point ordering so no need to store.
    forAll(patch_.localPoints(), pointI)
    {
        //label addedPointI =
        meshMod.addPoint
        (
            patch_.localPoints()[pointI],   // point
            pointToPointMap.size(),         // masterPointID
            -1,                             // zoneID
            true                            // inCell
        );
        pointToPointMap.append(pointI);

        //Pout<< "Added bottom point " << addedPointI
        //    << " at " << patch_.localPoints()[pointI]
        //    << "  from point " << pointI
        //    << endl;
    }


    // Introduce new points (one for every region)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList addedPoints(nLayers*regionPoints_.size());
    forAll(regionPoints_, regionI)
    {
        label pointI = regionPoints_[regionI];

        point pt = patch_.localPoints()[pointI];
        point disp = firstLayerDisp[regionI];
        for (label layerI = 0; layerI < nLayers; layerI++)
        {
            pt += disp;

            addedPoints[nLayers*regionI+layerI] = meshMod.addPoint
            (
                pt,                     // point
                pointToPointMap.size(), // masterPointID - used only addressing
                -1,                     // zoneID
                true                    // inCell
            );
            pointToPointMap.append(pointI);

            disp *= expansionRatio;
        }
    }


    // Add face on bottom side
    forAll(patch_.localFaces(), faceI)
    {
        meshMod.addFace
        (
            patch_.localFaces()[faceI].reverseFace(),// vertices
            addedCells[nLayers*faceI],  // own
            -1,                         // nei
            -1,                         // masterPointID
            -1,                         // masterEdgeID
            faceToFaceMap.size(),       // masterFaceID : current faceI
            true,                       // flipFaceFlux
            bottomPatchID[faceI],       // patchID
            -1,                         // zoneID
            false                       // zoneFlip
        );
        faceToFaceMap.append(-faceI-1); // points to flipped original face
        faceToEdgeMap.append(-1);

        //const face newF(patch_.localFaces()[faceI].reverseFace());
        //Pout<< "Added bottom face "
        //    << newF
        //    << " coords:" << UIndirectList<point>(meshMod.points(), newF)
        //    << " own " << addedCells[faceI]
        //    << " patch:" << bottomPatchID[faceI]
        //    << "  at " << patch_.faceCentres()[faceI]
        //    << endl;
    }

    // Add inbetween faces and face on top
    forAll(patch_.localFaces(), faceI)
    {
        // Get face in original ordering
        const face& f = patch_.localFaces()[faceI];

        face newF(f.size());

        for (label layerI = 0; layerI < nLayers; layerI++)
        {
            // Pick up point based on region and layer
            forAll(f, fp)
            {
                label region = pointRegions_[faceI][fp];
                newF[fp] = addedPoints[region*nLayers+layerI];
            }

            label own = addedCells[faceI*nLayers+layerI];
            label nei;
            label patchI;
            if (layerI == nLayers-1)
            {
                nei = -1;
                patchI = topPatchID[faceI];
            }
            else
            {
                nei = addedCells[faceI*nLayers+layerI+1];
                patchI = -1;
            }

            meshMod.addFace
            (
                newF,                       // vertices
                own,                        // own
                nei,                        // nei
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                faceToFaceMap.size(),       // masterFaceID : current faceI
                false,                      // flipFaceFlux
                patchI,                     // patchID
                -1,                         // zoneID
                false                       // zoneFlip
            );
            faceToFaceMap.append(faceI+1);  // unflipped
            faceToEdgeMap.append(-1);

            //Pout<< "Added inbetween face " << newF
            //    << " coords:" << UIndirectList<point>(meshMod.points(), newF)
            //    << " at layer " << layerI
            //    << " own " << own
            //    << " nei " << nei
            //    << "  at " << patch_.faceCentres()[faceI]
            //    << endl;
        }
    }


    // Add side faces
    // ~~~~~~~~~~~~~~
    // Note that we loop over edges multiple times so for edges with
    // two cyclic faces they get added in two passes (for correct ordering)

    // Pass1. Internal edges and first face of other edges
    forAll(extrudeEdgePatches, edgeI)
    {
        const labelList& eFaces = patch_.edgeFaces()[edgeI];
        const labelList& ePatches = extrudeEdgePatches[edgeI];

        if (ePatches.size() == 0)
        {
            // Internal face
            if (eFaces.size() != 2)
            {
                FatalErrorIn("createShellMesh::setRefinement(..)")
                    << "edge:" << edgeI
                    << " not internal but does not have side-patches defined."
                    << exit(FatalError);
            }
        }
        else
        {
            if (eFaces.size() != ePatches.size())
            {
                FatalErrorIn("createShellMesh::setRefinement(..)")
                    << "external/feature edge:" << edgeI
                    << " has " << eFaces.size() << " connected extruded faces "
                    << " but only " << ePatches.size()
                    << " boundary faces defined." << exit(FatalError);
            }
        }



        // Make face pointing in to eFaces[0] so out of new master face
        const face& f = patch_.localFaces()[eFaces[0]];
        const edge& e = patch_.edges()[edgeI];

        label fp0 = findIndex(f, e[0]);
        label fp1 = f.fcIndex(fp0);

        if (f[fp1] != e[1])
        {
            fp1 = fp0;
            fp0 = f.rcIndex(fp1);
        }

        face newF(4);

        for (label layerI = 0; layerI < nLayers; layerI++)
        {
            label region0 = pointRegions_[eFaces[0]][fp0];
            label region1 = pointRegions_[eFaces[0]][fp1];

            // Pick up points with correct normal
            if (layerI == 0)
            {
                newF[0] = f[fp0];
                newF[1] = f[fp1];
                newF[2] = addedPoints[nLayers*region1+layerI];
                newF[3] = addedPoints[nLayers*region0+layerI];
            }
            else
            {
                newF[0] = addedPoints[nLayers*region0+layerI-1];
                newF[1] = addedPoints[nLayers*region1+layerI-1];
                newF[2] = addedPoints[nLayers*region1+layerI];
                newF[3] = addedPoints[nLayers*region0+layerI];
            }

            // Optionally rotate so e[0] is always 0th vertex. Note that
            // this normally is automatically done by coupled face ordering
            // but with NOORDERING we have to do it ourselves.
            if (f[fp0] != e[0])
            {
                // rotate one back to get newF[1] (originating from e[0])
                // into newF[0]
                label v0 = newF[0];
                for (label i = 0; i < newF.size()-1; i++)
                {
                    newF[i] = newF[newF.fcIndex(i)];
                }
                newF.last() = v0;
            }


            label minCellI = addedCells[nLayers*eFaces[0]+layerI];
            label maxCellI;
            label patchI;
            if (ePatches.size() == 0)
            {
                maxCellI = addedCells[nLayers*eFaces[1]+layerI];
                if (minCellI > maxCellI)
                {
                    // Swap
                    Swap(minCellI, maxCellI);
                    newF = newF.reverseFace();
                }
                patchI = -1;
            }
            else
            {
                maxCellI = -1;
                patchI = ePatches[0];
            }

            //{
            //    Pout<< "Adding from face:" << patch_.faceCentres()[eFaces[0]]
            //        << " from edge:"
            //        << patch_.localPoints()[f[fp0]]
            //        << patch_.localPoints()[f[fp1]]
            //        << " at layer:" << layerI
            //        << " with new points:" << newF
            //        << " locations:"
            //        << UIndirectList<point>(meshMod.points(), newF)
            //        << " own:" << minCellI
            //        << " nei:" << maxCellI
            //        << endl;
            //}


            // newF already outwards pointing.
            meshMod.addFace
            (
                newF,                   // vertices
                minCellI,               // own
                maxCellI,               // nei
                -1,                     // masterPointID
                -1,                     // masterEdgeID
                faceToFaceMap.size(),   // masterFaceID
                false,                  // flipFaceFlux
                patchI,                 // patchID
                -1,                     // zoneID
                false                   // zoneFlip
            );
            faceToFaceMap.append(0);
            faceToEdgeMap.append(edgeI);
        }
    }

    // Pass2. Other faces of boundary edges
    forAll(extrudeEdgePatches, edgeI)
    {
        const labelList& eFaces = patch_.edgeFaces()[edgeI];
        const labelList& ePatches = extrudeEdgePatches[edgeI];

        if (ePatches.size() >= 2)
        {
            for (label i = 1; i < ePatches.size(); i++)
            {
                // Extrude eFaces[i]
                label minFaceI = eFaces[i];

                // Make face pointing in to eFaces[0] so out of new master face
                const face& f = patch_.localFaces()[minFaceI];

                const edge& e = patch_.edges()[edgeI];
                label fp0 = findIndex(f, e[0]);
                label fp1 = f.fcIndex(fp0);

                if (f[fp1] != e[1])
                {
                    fp1 = fp0;
                    fp0 = f.rcIndex(fp1);
                }

                face newF(4);
                for (label layerI = 0; layerI < nLayers; layerI++)
                {
                    label region0 = pointRegions_[minFaceI][fp0];
                    label region1 = pointRegions_[minFaceI][fp1];

                    if (layerI == 0)
                    {
                        newF[0] = f[fp0];
                        newF[1] = f[fp1];
                        newF[2] = addedPoints[nLayers*region1+layerI];
                        newF[3] = addedPoints[nLayers*region0+layerI];
                    }
                    else
                    {
                        newF[0] = addedPoints[nLayers*region0+layerI-1];
                        newF[1] = addedPoints[nLayers*region1+layerI-1];
                        newF[2] = addedPoints[nLayers*region1+layerI];
                        newF[3] = addedPoints[nLayers*region0+layerI];
                    }


                    // Optionally rotate so e[0] is always 0th vertex. Note that
                    // this normally is automatically done by coupled face
                    // ordering but with NOORDERING we have to do it ourselves.
                    if (f[fp0] != e[0])
                    {
                        // rotate one back to get newF[1] (originating
                        // from e[0]) into newF[0].
                        label v0 = newF[0];
                        for (label i = 0; i < newF.size()-1; i++)
                        {
                            newF[i] = newF[newF.fcIndex(i)];
                        }
                        newF.last() = v0;
                    }
                    ////if (ePatches.size() == 0)
                    //{
                    //    Pout<< "Adding from MULTI face:"
                    //        << patch_.faceCentres()[minFaceI]
                    //        << " from edge:"
                    //        << patch_.localPoints()[f[fp0]]
                    //        << patch_.localPoints()[f[fp1]]
                    //        << " at layer:" << layerI
                    //        << " with new points:" << newF
                    //        << " locations:"
                    //        << UIndirectList<point>(meshMod.points(), newF)
                    //        << endl;
                    //}

                    // newF already outwards pointing.
                    meshMod.addFace
                    (
                        newF,                   // vertices
                        addedCells[nLayers*minFaceI+layerI],   // own
                        -1,                     // nei
                        -1,                     // masterPointID
                        -1,                     // masterEdgeID
                        faceToFaceMap.size(),   // masterFaceID
                        false,                  // flipFaceFlux
                        ePatches[i],            // patchID
                        -1,                     // zoneID
                        false                   // zoneFlip
                    );
                    faceToFaceMap.append(0);
                    faceToEdgeMap.append(edgeI);
                }
            }
        }
    }


    cellToFaceMap_.transfer(cellToFaceMap);
    faceToFaceMap_.transfer(faceToFaceMap);
    faceToEdgeMap_.transfer(faceToEdgeMap);
    pointToPointMap_.transfer(pointToPointMap);
}


void Foam::createShellMesh::updateMesh(const mapPolyMesh& map)
{
    inplaceReorder(map.reverseCellMap(), cellToFaceMap_);
    inplaceReorder(map.reverseFaceMap(), faceToFaceMap_);
    inplaceReorder(map.reverseFaceMap(), faceToEdgeMap_);
    inplaceReorder(map.reversePointMap(), pointToPointMap_);
}


// ************************************************************************* //

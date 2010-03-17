/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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
    const pointField& thickness,
    const labelList& extrudeMasterPatchID,
    const labelList& extrudeSlavePatchID,
    const labelListList& extrudeEdgePatches,
    polyTopoChange& meshMod
)
{
    if (thickness.size() != regionPoints_.size())
    {
        FatalErrorIn("createShellMesh::setRefinement(..)")
            << "nRegions:" << regionPoints_.size()
            << " thickness:" << thickness.size()
            << exit(FatalError);
    }

    if
    (
        extrudeMasterPatchID.size() != patch_.size()
     && extrudeSlavePatchID.size() != patch_.size()
    )
    {
        FatalErrorIn("createShellMesh::setRefinement(..)")
            << "nFaces:" << patch_.size()
            << " extrudeMasterPatchID:" << extrudeMasterPatchID.size()
            << " extrudeSlavePatchID:" << extrudeSlavePatchID.size()
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
    DynamicList<label> cellToFaceMap(patch_.size());
    // From face to patch+turning index
    DynamicList<label> faceToFaceMap(2*patch_.size()+patch_.nEdges());
    // From face to patch edge index
    DynamicList<label> faceToEdgeMap(patch_.nEdges()+patch_.nEdges());
    // From point to patch point index
    DynamicList<label> pointToPointMap(2*patch_.nPoints());


    // Introduce new cell for every face
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList addedCells(patch_.size());
    forAll(patch_, faceI)
    {
        addedCells[faceI] = meshMod.addCell
        (
            -1,                     // masterPointID
            -1,                     // masterEdgeID
            -1,                     // masterFaceID
            cellToFaceMap.size(),   // masterCellID
            -1                      // zoneID
        );
        cellToFaceMap.append(faceI);
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

    labelList addedPoints(regionPoints_.size());
    forAll(regionPoints_, regionI)
    {
        label pointI = regionPoints_[regionI];
        point extrudedPt = patch_.localPoints()[pointI] + thickness[regionI];

        addedPoints[regionI] = meshMod.addPoint
        (
            extrudedPt,             // point
            pointToPointMap.size(), // masterPointID - used only addressing
            -1,                     // zoneID
            true                    // inCell
        );
        pointToPointMap.append(pointI);

        //Pout<< "Added top point " << addedPoints[regionI]
        //    << " at " << extrudedPt
        //    << "  from point " << pointI
        //    << endl;
    }


    // Add face on patch' master side
    //labelList masterFaces(patch_.localFaces().size());
    forAll(patch_.localFaces(), faceI)
    {
        //masterFaces[faceI] =
        meshMod.addFace
        (
            patch_.localFaces()[faceI].reverseFace(),// vertices
            addedCells[faceI],          // own
            -1,                         // nei
            -1,                         // masterPointID
            -1,                         // masterEdgeID
            faceToFaceMap.size(),       // masterFaceID : current faceI
            true,                       // flipFaceFlux
            extrudeMasterPatchID[faceI],// patchID
            -1,                         // zoneID
            false                       // zoneFlip
        );
        faceToFaceMap.append(faceI+1);  // points to unflipped original face
        faceToEdgeMap.append(-1);

        //Pout<< "Added master face "
        //    << patch_.localFaces()[faceI].reverseFace()
        //    << " own " << addedCells[faceI]
        //    << "  at " << patch_.faceCentres()[faceI]
        //    << endl;

    }

    // Add face on patch' slave side
    //labelList slaveFaces(patch_.localFaces().size());
    forAll(patch_.localFaces(), faceI)
    {
        // Get face in original ordering
        const face& f = patch_.localFaces()[faceI];

        // Pick up point based on region
        face newF(f.size());
        forAll(f, fp)
        {
            label region = pointRegions_[faceI][fp];
            newF[fp] = addedPoints[region];
        }

        //slaveFaces[faceI] =
        meshMod.addFace
        (
            newF,                       // vertices
            addedCells[faceI],          // own
            -1,                         // nei
            -1,                         // masterPointID
            -1,                         // masterEdgeID
            faceToFaceMap.size(),       // masterFaceID : current faceI
            false,                      // flipFaceFlux
            extrudeSlavePatchID[faceI], // patchID
            -1,                         // zoneID
            false                       // zoneFlip
        );
        faceToFaceMap.append(-faceI-1);
        faceToEdgeMap.append(-1);

        //Pout<< "Added slave face " << newF
        //    << " own " << addedCells[faceI]
        //    << "  at " << patch_.faceCentres()[faceI]
        //    << endl;
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
            // internal face.

            if (eFaces.size() != 2)
            {
                FatalErrorIn("createShellMesh::setRefinement(..)")
                    << "edge:" << edgeI
                    << " not internal but does not have side-patches defined."
                    << exit(FatalError);
            }

            // Extrude

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
            newF[0] = f[fp0];
            newF[1] = f[fp1];
            newF[2] = addedPoints[pointRegions_[eFaces[0]][fp1]];
            newF[3] = addedPoints[pointRegions_[eFaces[0]][fp0]];

            label minCellI = addedCells[eFaces[0]];
            label maxCellI = addedCells[eFaces[1]];

            if (minCellI > maxCellI)
            {
                // Swap
                Swap(minCellI, maxCellI);
                newF = newF.reverseFace();
            }

            //Pout<< "for internal edge:" << e
            //    << " at:" << patch_.localPoints()[e[0]]
            //    << patch_.localPoints()[e[1]]
            //    << " adding face:" << newF
            //    << " from f:" << f
            //    << " inbetween " << minCellI << " and " << maxCellI << endl;

            // newF already outwards pointing.
            meshMod.addFace
            (
                newF,       // vertices
                minCellI,   // own
                maxCellI,   // nei
                -1,         // masterPointID
                -1,         // masterEdgeID
                faceToFaceMap.size(),   // masterFaceID
                false,      // flipFaceFlux
                -1,         // patchID
                -1,         // zoneID
                false       // zoneFlip
            );
            faceToFaceMap.append(0);
            faceToEdgeMap.append(edgeI);
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


            // Extrude eFaces[0]
            label minFaceI = eFaces[0];

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
            newF[0] = f[fp0];
            newF[1] = f[fp1];
            newF[2] = addedPoints[pointRegions_[minFaceI][fp1]];
            newF[3] = addedPoints[pointRegions_[minFaceI][fp0]];


            //Pout<< "for external edge:" << e
            //    << " at:" << patch_.localPoints()[e[0]]
            //    << patch_.localPoints()[e[1]]
            //    << " adding first patch face:" << newF
            //    << " from:" << f
            //    << " into patch:" << ePatches[0]
            //    << " own:" << addedCells[minFaceI]
            //    << endl;


            // newF already outwards pointing.
            meshMod.addFace
            (
                newF,                   // vertices
                addedCells[minFaceI],   // own
                -1,                     // nei
                -1,                     // masterPointID
                -1,                     // masterEdgeID
                faceToFaceMap.size(),   // masterFaceID
                false,                  // flipFaceFlux
                ePatches[0],            // patchID
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
                newF[0] = f[fp0];
                newF[1] = f[fp1];
                newF[2] = addedPoints[pointRegions_[minFaceI][fp1]];
                newF[3] = addedPoints[pointRegions_[minFaceI][fp0]];

                //Pout<< "for external edge:" << e
                //    << " at:" << patch_.localPoints()[e[0]]
                //    << patch_.localPoints()[e[1]]
                //    << " adding patch face:" << newF
                //    << " from:" << f
                //    << " into patch:" << ePatches[i]
                //    << endl;

                // newF already outwards pointing.
                meshMod.addFace
                (
                    newF,                   // vertices
                    addedCells[minFaceI],   // own
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

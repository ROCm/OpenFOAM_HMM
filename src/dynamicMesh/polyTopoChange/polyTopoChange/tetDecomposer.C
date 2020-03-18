/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "tetDecomposer.H"
#include "meshTools.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "OFstream.H"
#include "edgeHashes.H"
#include "syncTools.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tetDecomposer, 0);
}

const Foam::Enum
<
    Foam::tetDecomposer::decompositionType
>
Foam::tetDecomposer::decompositionTypeNames
({
    { decompositionType::FACE_CENTRE_TRIS,  "faceCentre" },
    { decompositionType::FACE_DIAG_TRIS, "faceDiagonal" },
    { decompositionType::PYRAMID, "pyramid" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::tetDecomposer::modifyFace
(
    polyTopoChange& meshMod,
    const face& f,
    const label facei,
    const label own,
    const label nei,
    const label patchi,
    const label zoneI,
    const bool zoneFlip
) const
{
    // First usage of face. Modify.
    if (nei == -1 || own < nei)
    {
        meshMod.modifyFace
        (
            f,                          // modified face
            facei,                      // label of face
            own,                        // owner
            nei,                        // neighbour
            false,                      // face flip
            patchi,                     // patch for face
            zoneI,                      // zone for face
            zoneFlip                    // face flip in zone
        );
    }
    else
    {
        meshMod.modifyFace
        (
            f.reverseFace(),            // modified face
            facei,                      // label of face
            nei,                        // owner
            own,                        // neighbour
            true,                       // face flip
            patchi,                     // patch for face
            zoneI,                      // zone for face
            !zoneFlip                   // face flip in zone
        );
    }
}


void Foam::tetDecomposer::addFace
(
    polyTopoChange& meshMod,
    const face& f,
    const label own,
    const label nei,
    const label masterPointID,
    const label masterEdgeID,
    const label masterFaceID,
    const label patchi,
    const label zoneI,
    const bool zoneFlip
) const
{
    // Second or more usage of face. Add.
    if (nei == -1 || own < nei)
    {
        meshMod.addFace
        (
            f,                          // modified face
            own,                        // owner
            nei,                        // neighbour
            masterPointID,              // master point
            masterEdgeID,               // master edge
            masterFaceID,               // master face
            false,                      // face flip
            patchi,                     // patch for face
            zoneI,                      // zone for face
            zoneFlip                    // face flip in zone
        );
    }
    else
    {
        meshMod.addFace
        (
            f.reverseFace(),            // modified face
            nei,                        // owner
            own,                        // neighbour
            masterPointID,              // master point
            masterEdgeID,               // master edge
            masterFaceID,               // master face
            true,                       // face flip
            patchi,                     // patch for face
            zoneI,                      // zone for face
            !zoneFlip                   // face flip in zone
        );
    }
}


// Work out triangle index given the starting vertex in the face
Foam::label Foam::tetDecomposer::triIndex(const label facei, const label fp)
const
{
    const face& f = mesh_.faces()[facei];
    const label fp0 = max(0, mesh_.tetBasePtIs()[facei]);

    // Work out triangle index on this face
    label thisTrii;
    if (fp == fp0)
    {
        thisTrii = 0;
    }
    else if (fp == f.rcIndex(fp0))
    {
        thisTrii = f.size()-3;
    }
    else
    {
        thisTrii = (fp-fp0-1) % (f.size()-2);
    }
    return thisTrii;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetDecomposer::tetDecomposer(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tetDecomposer::setRefinement
(
    const decompositionType decomposeType,
    const bitSet& decomposeCell,
    polyTopoChange& meshMod
)
{
    cellToPoint_.setSize(mesh_.nCells(), -1);
    forAll(mesh_.cellCentres(), celli)
    {
        if (decomposeCell[celli])
        {
            // Any point on the cell
            label masterPointi = mesh_.faces()[mesh_.cells()[celli][0]][0];

            cellToPoint_[celli] = meshMod.addPoint
            (
                mesh_.cellCentres()[celli],
                masterPointi,
                -1,
                true
            );
        }
    }


    // Determine for every face whether it borders a cell that is decomposed
    bitSet decomposeFace(mesh_.nFaces());
    {
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label nei = mesh_.faceNeighbour()[facei];
            if (decomposeCell[own] || decomposeCell[nei])
            {
                decomposeFace[facei] = true;
            }
        }

        boolList neiDecomposeCell(mesh_.nBoundaryFaces());
        forAll(neiDecomposeCell, bFacei)
        {
            label facei = mesh_.nInternalFaces()+bFacei;
            label own = mesh_.faceOwner()[facei];
            neiDecomposeCell[bFacei] = decomposeCell[own];
        }
        syncTools::swapBoundaryFaceList(mesh_, neiDecomposeCell);

        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            label own = mesh_.faceOwner()[facei];
            label bFacei = facei-mesh_.nInternalFaces();
            if (decomposeCell[own] || neiDecomposeCell[bFacei])
            {
                decomposeFace[facei] = true;
            }
        }
    }


    // Add face centre points
    if (decomposeType == FACE_CENTRE_TRIS)
    {
        faceToPoint_.setSize(mesh_.nFaces(), -1);
        forAll(mesh_.faceCentres(), facei)
        {
            if (decomposeFace[facei])
            {
                // Any point on the face
                const label masterPointi = mesh_.faces()[facei][0];

                faceToPoint_[facei] = meshMod.addPoint
                (
                    mesh_.faceCentres()[facei],
                    masterPointi,
                    -1,
                    true
                );
            }
        }
    }


    // Per face, per point (faceCentre) or triangle (faceDiag) the (existing
    // or added) cell on either side
    faceOwnerCells_.setSize(mesh_.nFaces());
    faceNeighbourCells_.setSize(mesh_.nFaces());

    if (decomposeType == FACE_CENTRE_TRIS)
    {
        forAll(faceOwnerCells_, facei)
        {
            const face& f = mesh_.faces()[facei];
            faceOwnerCells_[facei].setSize(f.size(), -1);
            faceNeighbourCells_[facei].setSize(f.size(), -1);
        }
    }
    else if (decomposeType == FACE_DIAG_TRIS)
    {
        // Force construction of diagonal decomposition
        (void)mesh_.tetBasePtIs();

        forAll(faceOwnerCells_, facei)
        {
            const face& f = mesh_.faces()[facei];
            faceOwnerCells_[facei].setSize(f.size()-2, -1);
            faceNeighbourCells_[facei].setSize(f.size()-2, -1);
        }
    }
    else
    {
        forAll(faceOwnerCells_, facei)
        {
            faceOwnerCells_[facei].setSize(1, -1);
            faceNeighbourCells_[facei].setSize(1, -1);
        }
    }


    // Add internal cells. Note: done in same order as pyramid triangle
    // creation later to maintain same ordering.
    forAll(mesh_.cells(), celli)
    {
        const cell& cFaces = mesh_.cells()[celli];

        // Whether cell has already been modified (all other cells get added)
        bool modifiedCell = false;

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];
            const face& f = mesh_.faces()[facei];

            // Get reference to either owner or neighbour
            labelList& added =
            (
                (mesh_.faceOwner()[facei] == celli)
              ? faceOwnerCells_[facei]
              : faceNeighbourCells_[facei]
            );

            if (decomposeCell[celli])
            {
                if (decomposeType == FACE_CENTRE_TRIS)
                {
                    forAll(f, fp)
                    {
                        if (!modifiedCell)
                        {
                            // Reuse cell itself
                            added[fp] = celli;
                            modifiedCell = true;
                        }
                        else
                        {
                            added[fp] = meshMod.addCell
                            (
                                -1,     // masterPoint
                                -1,     // masterEdge
                                -1,     // masterFace
                                celli,  // masterCell
                                mesh_.cellZones().whichZone(celli)
                            );
                        }
                    }
                }
                else if (decomposeType == FACE_DIAG_TRIS)
                {
                    for (label triI = 0; triI < f.size()-2; triI++)
                    {
                        if (!modifiedCell)
                        {
                            // Reuse cell itself
                            added[triI] = celli;
                            modifiedCell = true;
                        }
                        else
                        {
                            added[triI] = meshMod.addCell
                            (
                                -1,     // masterPoint
                                -1,     // masterEdge
                                -1,     // masterFace
                                celli,  // masterCell
                                mesh_.cellZones().whichZone(celli)
                            );
                            //Pout<< "For cell:" << celli
                            //    << " at:" << mesh_.cellCentres()[celli]
                            //    << " face:" << facei
                            //    << " at:" << mesh_.faceCentres()[facei]
                            //    << " tri:" << triI
                            //    << " added cell:" << added[triI] << endl;
                        }
                    }
                }
                else // if (decomposeType == PYRAMID)
                {
                    // Pyramidal decomposition.
                    // Assign same cell to all slots
                    if (!modifiedCell)
                    {
                        // Reuse cell itself
                        added = celli;
                        modifiedCell = true;
                    }
                    else
                    {
                        added = meshMod.addCell
                        (
                            -1,     // masterPoint
                            -1,     // masterEdge
                            -1,     // masterFace
                            celli,  // masterCell
                            mesh_.cellZones().whichZone(celli)
                        );
                    }
                }
            }
            else
            {
                // All vertices/triangles address to original cell
                added = celli;
            }
        }
    }



    // Add triangle faces
    face triangle(3);

    forAll(mesh_.faces(), facei)
    {
        label own = mesh_.faceOwner()[facei];
        const labelList& addedOwn = faceOwnerCells_[facei];
        const labelList& addedNei = faceNeighbourCells_[facei];
        const face& f = mesh_.faces()[facei];

        label patchi = -1;
        if (facei >= mesh_.nInternalFaces())
        {
            patchi = mesh_.boundaryMesh().whichPatch(facei);
        }

        label zoneI = mesh_.faceZones().whichZone(facei);
        bool zoneFlip = false;
        if (zoneI != -1)
        {
            const faceZone& fz = mesh_.faceZones()[zoneI];
            zoneFlip = fz.flipMap()[fz.whichFace(facei)];
        }

        //Pout<< "Face:" << facei << " at:" << mesh_.faceCentres()[facei]
        //    << endl;

        if (decomposeType == FACE_CENTRE_TRIS && decomposeFace[facei])
        {
            forAll(f, fp)
            {
                // 1. Front triangle (decomposition of face itself)
                //    (between owner and neighbour cell)
                {
                    triangle[0] = f[fp];
                    triangle[1] = f[f.fcIndex(fp)];
                    triangle[2] = faceToPoint_[facei];

                    //Pout<< "    triangle:" << triangle
                    //    << " points:"
                    //    << UIndirectList<point>(meshMod.points(), triangle)
                    //    << " between:" << addedOwn[fp]
                    //    << " and:" << addedNei[fp] << endl;


                    if (fp == 0)
                    {
                        modifyFace
                        (
                            meshMod,
                            triangle,
                            facei,
                            addedOwn[fp],
                            addedNei[fp],
                            patchi,
                            zoneI,
                            zoneFlip
                        );
                    }
                    else
                    {
                        addFace
                        (
                            meshMod,
                            triangle,
                            addedOwn[fp],
                            addedNei[fp],
                            -1,                 //point
                            -1,                 //edge
                            facei,              //face
                            patchi,
                            zoneI,
                            zoneFlip
                        );
                    }
                }


                // 2. Within owner cell - to cell centre
                if (decomposeCell[own])
                {
                    label newOwn = addedOwn[f.rcIndex(fp)];
                    label newNei = addedOwn[fp];

                    triangle[0] = f[fp];
                    triangle[1] = cellToPoint_[own];
                    triangle[2] = faceToPoint_[facei];

                    addFace
                    (
                        meshMod,
                        triangle,
                        newOwn,
                        newNei,
                        f[fp],      //point
                        -1,         //edge
                        -1,         //face
                        -1,         //patchi
                        -1,         //zone
                        false
                    );
                }
                // 2b. Within neighbour cell - to cell centre
                if
                (
                    facei < mesh_.nInternalFaces()
                 && decomposeCell[mesh_.faceNeighbour()[facei]]
                )
                {
                    label newOwn = addedNei[f.rcIndex(fp)];
                    label newNei = addedNei[fp];

                    triangle[0] = f[fp];
                    triangle[1] = faceToPoint_[facei];
                    triangle[2] = cellToPoint_[mesh_.faceNeighbour()[facei]];

                    addFace
                    (
                        meshMod,
                        triangle,
                        newOwn,
                        newNei,
                        f[fp],      //point
                        -1,         //edge
                        -1,         //face
                        -1,         //patchi
                        -1,         //zone
                        false
                    );
                }
            }
        }
        else if (decomposeType == FACE_DIAG_TRIS && decomposeFace[facei])
        {
            label fp0 = max(mesh_.tetBasePtIs()[facei], 0);
            label fp = f.fcIndex(fp0);

            for (label triI = 0; triI < f.size()-2; triI++)
            {
                label nextTri = triI+1;
                if (nextTri >= f.size()-2)
                {
                    nextTri -= f.size()-2;
                }
                label nextFp = f.fcIndex(fp);


                // Triangle triI consisiting of f[fp0], f[fp], f[nextFp]


                // 1. Front triangle (decomposition of face itself)
                //    (between owner and neighbour cell)
                {
                    triangle[0] = f[fp0];
                    triangle[1] = f[fp];
                    triangle[2] = f[nextFp];

                    if (triI == 0)
                    {
                        modifyFace
                        (
                            meshMod,
                            triangle,
                            facei,
                            addedOwn[triI],
                            addedNei[triI],
                            patchi,
                            zoneI,
                            zoneFlip
                        );
                    }
                    else
                    {
                        addFace
                        (
                            meshMod,
                            triangle,
                            addedOwn[triI],
                            addedNei[triI],
                            -1,                 //point
                            -1,                 //edge
                            facei,              //face
                            patchi,
                            zoneI,
                            zoneFlip
                        );
                    }
                }


                // 2. Within owner cell - diagonal to cell centre
                if (triI < f.size()-3)
                {
                    if (decomposeCell[own])
                    {
                        label newOwn = addedOwn[triI];
                        label newNei = addedOwn[nextTri];

                        triangle[0] = f[fp0];
                        triangle[1] = f[nextFp];
                        triangle[2] = cellToPoint_[own];

                        addFace
                        (
                            meshMod,
                            triangle,
                            newOwn,
                            newNei,
                            f[fp],      //point
                            -1,         //edge
                            -1,         //face
                            -1,         //patchi
                            -1,         //zone
                            false
                        );
                    }

                    // 2b. Within neighbour cell - to cell centre
                    if
                    (
                        facei < mesh_.nInternalFaces()
                     && decomposeCell[mesh_.faceNeighbour()[facei]]
                    )
                    {
                        label newOwn = addedNei[triI];
                        label newNei = addedNei[nextTri];

                        triangle[0] = f[nextFp];
                        triangle[1] = f[fp0];
                        triangle[2] =
                            cellToPoint_[mesh_.faceNeighbour()[facei]];

                        addFace
                        (
                            meshMod,
                            triangle,
                            newOwn,
                            newNei,
                            f[fp],      //point
                            -1,         //edge
                            -1,         //face
                            -1,         //patchi
                            -1,         //zone
                            false
                        );
                    }
                }


                fp = nextFp;
            }
        }
        else
        {
            // No decomposition. Use zero'th slot.
            modifyFace
            (
                meshMod,
                f,
                facei,
                addedOwn[0],
                addedNei[0],
                patchi,
                zoneI,
                zoneFlip
            );
        }
    }



    // Add triangles for all edges.
    EdgeMap<label> edgeToFace;

    forAll(mesh_.cells(), celli)
    {
        const cell& cFaces = mesh_.cells()[celli];

        edgeToFace.clear();

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];

            if (decomposeCell[celli])
            {
                const face& f = mesh_.faces()[facei];
                //const labelList& fEdges = mesh_.faceEdges()[facei];
                forAll(f, fp)
                {
                    label p0 = f[fp];
                    label p1 = f[f.fcIndex(fp)];
                    const edge e(p0, p1);

                    EdgeMap<label>::const_iterator edgeFnd = edgeToFace.find(e);
                    if (edgeFnd == edgeToFace.end())
                    {
                        edgeToFace.insert(e, facei);
                    }
                    else
                    {
                        // Found the other face on the edge.
                        label otherFacei = edgeFnd();
                        const face& otherF = mesh_.faces()[otherFacei];

                        // Found the other face on the edge. Note that since
                        // we are looping in the same order the tets added for
                        // otherFacei will be before those of facei

                        label otherFp = otherF.find(p0);
                        if (otherF.nextLabel(otherFp) == p1)
                        {
                            // ok. otherFp is first vertex of edge.
                        }
                        else if (otherF.prevLabel(otherFp) == p1)
                        {
                            otherFp = otherF.rcIndex(otherFp);
                        }
                        else
                        {
                            FatalErrorInFunction
                                << "problem." << abort(FatalError);
                        }


                        // Triangle from edge to cell centre
                        if (mesh_.faceOwner()[facei] == celli)
                        {
                            triangle[0] = p0;
                            triangle[1] = p1;
                            triangle[2] = cellToPoint_[celli];
                        }
                        else
                        {
                            triangle[0] = p1;
                            triangle[1] = p0;
                            triangle[2] = cellToPoint_[celli];
                        }

                        // Determine tets on either side
                        label thisTet, otherTet;

                        if (decomposeType == FACE_CENTRE_TRIS)
                        {
                            if (mesh_.faceOwner()[facei] == celli)
                            {
                                thisTet = faceOwnerCells_[facei][fp];
                            }
                            else
                            {
                                thisTet = faceNeighbourCells_[facei][fp];
                            }

                            if (mesh_.faceOwner()[otherFacei] == celli)
                            {
                                otherTet =
                                    faceOwnerCells_[otherFacei][otherFp];
                            }
                            else
                            {
                                otherTet =
                                    faceNeighbourCells_[otherFacei][otherFp];
                            }
                        }
                        else if (decomposeType == FACE_DIAG_TRIS)
                        {
                            label thisTriI = triIndex(facei, fp);
                            if (mesh_.faceOwner()[facei] == celli)
                            {
                                thisTet = faceOwnerCells_[facei][thisTriI];
                            }
                            else
                            {
                                thisTet = faceNeighbourCells_[facei][thisTriI];
                            }

                            label otherTriI = triIndex(otherFacei, otherFp);
                            if (mesh_.faceOwner()[otherFacei] == celli)
                            {
                                otherTet =
                                    faceOwnerCells_[otherFacei][otherTriI];
                            }
                            else
                            {
                                otherTet =
                                    faceNeighbourCells_[otherFacei][otherTriI];
                            }
                        }
                        else
                        {
                            if (mesh_.faceOwner()[facei] == celli)
                            {
                                thisTet = faceOwnerCells_[facei][0];
                            }
                            else
                            {
                                thisTet = faceNeighbourCells_[facei][0];
                            }
                            if (mesh_.faceOwner()[otherFacei] == celli)
                            {
                                otherTet = faceOwnerCells_[otherFacei][0];
                            }
                            else
                            {
                                otherTet =
                                    faceNeighbourCells_[otherFacei][0];
                            }
                        }

                        addFace
                        (
                            meshMod,
                            triangle,
                            otherTet,
                            thisTet,
                            -1,         //masterPoint
                            -1,         //fEdges[fp], //masterEdge
                            facei,      //masterFace
                            -1,         //patchi
                            -1,         //zone
                            false
                        );
                    }
                }
            }
        }
    }
}


void Foam::tetDecomposer::updateMesh(const mapPolyMesh& map)
{
    inplaceRenumber(map.reversePointMap(), cellToPoint_);
    inplaceRenumber(map.reversePointMap(), faceToPoint_);

    forAll(faceOwnerCells_, facei)
    {
        inplaceRenumber(map.reverseCellMap(), faceOwnerCells_[facei]);
    }
    forAll(faceNeighbourCells_, facei)
    {
        inplaceRenumber(map.reverseCellMap(), faceNeighbourCells_[facei]);
    }
}


// ************************************************************************* //

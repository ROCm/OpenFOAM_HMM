/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "CV3D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CV3D::calcDualMesh
(
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts
)
{
    // ~~~~~~~~~~~ removing short edges by indexing dual vertices ~~~~~~~~~~~~~~

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        cit->cellIndex() = -1;
    }

    points.setSize(number_of_cells());

    label dualVerti = 0;

    // Scanning by number of short (dual) edges (nSE) attached to the
    // circumcentre of each Delaunay tet.  A Delaunay tet may only have four
    // dual edges emanating from its circumcentre, assigning positions and
    // indices to those with 4 short edges attached first, then >= 3, then >= 2
    // etc.
    for (label nSE = 4; nSE >= 0; nSE--)
    {
        Info<< nl << "Scanning for dual vertices with >= "
            << nSE
            << " short edges attached." << endl;

        for
        (
            Triangulation::Finite_cells_iterator cit = finite_cells_begin();
            cit != finite_cells_end();
            ++cit
        )
        {
            // If the Delaunay tet has an index already then it has either
            // evaluated itself and taken action or has had its index dictated
            // by a neighbouring tet with more short edges attached.

            if (cit->cellIndex() == -1)
            {
                point dualVertex = topoint(dual(cit));

                label shortEdges = 0;

                List<bool> edgeIsShort(4, false);

                List<bool> neighbourAlreadyIndexed(4, false);

                // Loop over the four facets of the Delaunay tet
                for (label f = 0; f < 4; f++)
                {
                    // Check that at least one of the vertices of the facet is
                    // an internal or boundary point
                    if
                    (
                        cit->vertex(vertex_triple_index(f, 0))->
                        internalOrBoundaryPoint()
                     || cit->vertex(vertex_triple_index(f, 1))->
                        internalOrBoundaryPoint()
                     || cit->vertex(vertex_triple_index(f, 2))->
                        internalOrBoundaryPoint()
                    )
                    {
                        point neighDualVertex;

                        label cNI = cit->neighbor(f)->cellIndex();

                        if (cNI == -1)
                        {
                            neighDualVertex = topoint(dual(cit->neighbor(f)));
                        }
                        else
                        {
                            neighDualVertex = points[cNI];
                        }

                        if
                        (
                            magSqr(dualVertex - neighDualVertex)
                            < tols_.minEdgeLen2
                        )
                        {
                            edgeIsShort[f] = true;

                            if (cNI > -1)
                            {
                                neighbourAlreadyIndexed[f] = true;
                            }

                            shortEdges++;
                        }
                    }
                }

                if (nSE == 0 && shortEdges == 0)
                {
                    // Final iteration and no short edges are found, index
                    // remaining dual vertices.

                    if
                    (
                        cit->vertex(0)->internalOrBoundaryPoint()
                     || cit->vertex(1)->internalOrBoundaryPoint()
                     || cit->vertex(2)->internalOrBoundaryPoint()
                     || cit->vertex(3)->internalOrBoundaryPoint()
                    )
                    {
                        cit->cellIndex() = dualVerti;
                        points[dualVerti] = dualVertex;
                        dualVerti++;
                    }
                }
                else if
                (
                    shortEdges >= nSE
                )
                {
                    // Info<< neighbourAlreadyIndexed << ' '
                    //     << edgeIsShort << endl;

                    label numUnindexedNeighbours = 1;

                    for (label f = 0; f < 4; f++)
                    {
                        if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
                        {
                            dualVertex += topoint(dual(cit->neighbor(f)));

                            numUnindexedNeighbours++;
                        }
                    }

                    dualVertex /= numUnindexedNeighbours;

                    label nearestExistingIndex = -1;

                    point nearestIndexedNeighbourPos = vector::zero;

                    scalar minDistSqrToNearestIndexedNeighbour = VGREAT;

                    for (label f = 0; f < 4; f++)
                    {
                        if (edgeIsShort[f] && neighbourAlreadyIndexed[f])
                        {
                            label cNI = cit->neighbor(f)->cellIndex();

                            point indexedNeighbourPos = points[cNI];

                            if
                            (
                                magSqr(indexedNeighbourPos - dualVertex)
                                < minDistSqrToNearestIndexedNeighbour
                            )
                            {
                                nearestExistingIndex = cNI;

                                nearestIndexedNeighbourPos =
                                indexedNeighbourPos;

                                minDistSqrToNearestIndexedNeighbour =
                                magSqr(indexedNeighbourPos - dualVertex);
                            }
                        }
                    }

                    if
                    (
                        nearestExistingIndex > -1
                     && minDistSqrToNearestIndexedNeighbour < tols_.minEdgeLen2
                    )
                    {
                        points[nearestExistingIndex] =
                        0.5*(dualVertex + nearestIndexedNeighbourPos);

                        for (label f = 0; f < 4; f++)
                        {
                            if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
                            {
                                cit->neighbor(f)->cellIndex() =
                                nearestExistingIndex;
                            }
                        }

                        cit->cellIndex() = nearestExistingIndex;
                    }
                    else
                    {
                        for (label f = 0; f < 4; f++)
                        {
                            if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
                            {
                                cit->neighbor(f)->cellIndex() = dualVerti;
                            }
                        }

                        cit->cellIndex() = dualVerti;

                        points[dualVerti] = dualVertex;

                        dualVerti++;
                    }
                }
            }
        }
    }

    points.setSize(dualVerti);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~ dual cell indexing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // assigns an index to the Delaunay vertices which will be the dual cell
    // index used for owner neighbour assignment.

    // The indices of the points are reset which destroys the point-pair
    // matching, so the type of each vertex are reset to avoid any ambiguity.

    label dualCelli = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            vit->type() = Vb::INTERNAL_POINT;
            vit->index() = dualCelli;
            dualCelli++;
        }
        else
        {
            vit->type() = Vb::FAR_POINT;
            vit->index() = -1;
        }
    }

    // ~~~~~~~~~~~~ dual face and owner neighbour construction ~~~~~~~~~~~~~~~~~

    label nPatches = qSurf_.patches().size() + 1;

    label defaultPatchIndex = qSurf_.patches().size();

    patchNames.setSize(nPatches);

    const geometricSurfacePatchList& surfacePatches = qSurf_.patches();

    forAll(surfacePatches, sP)
    {
        patchNames[sP] = surfacePatches[sP].name();
    }

    patchNames[defaultPatchIndex] = "CV3D_default_patch";

    patchSizes.setSize(nPatches);

    patchStarts.setSize(nPatches);

    List<DynamicList<face> > patchFaces(nPatches, DynamicList<face>(0));

    List<DynamicList<label> > patchOwners(nPatches, DynamicList<label>(0));

    faces.setSize(number_of_edges());

    owner.setSize(number_of_edges());

    neighbour.setSize(number_of_edges());

    label dualFacei = 0;

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Cell_handle c = eit->first;
        Vertex_handle vA = c->vertex(eit->second);
        Vertex_handle vB = c->vertex(eit->third);

        if
        (
            vA->internalOrBoundaryPoint()
         || vB->internalOrBoundaryPoint()
        )
        {
            Cell_circulator ccStart = incident_cells(*eit);
            Cell_circulator cc1 = ccStart;
            Cell_circulator cc2 = cc1;

            // Advance the second circulator so that it always stays on the next
            // cell around the edge;
            cc2++;

            DynamicList<label> verticesOnFace;

            do
            {
                label cc1I = cc1->cellIndex();

                label cc2I = cc2->cellIndex();


                if (cc1I < 0 || cc2I < 0)
                {
                    FatalErrorIn("Foam::CV3D::calcDualMesh")
                        << "Dual face uses circumcenter defined by a "
                        << "Delaunay tetrahedron with no internal "
                        << "or boundary points.  Defining Delaunay edge ends: "
                        << topoint(vA->point()) << " "
                        << topoint(vB->point()) << nl
                        << exit(FatalError);
                }

                if (cc1I != cc2I)
                {
                    verticesOnFace.append(cc1I);
                }

                cc1++;

                cc2++;
            } while (cc1 != ccStart);

            verticesOnFace.shrink();

            if (verticesOnFace.size() >= 3)
            {
                face newDualFace(verticesOnFace);

                label dcA = vA->index();

                if (!vA->internalOrBoundaryPoint())
                {
                    dcA = -1;
                }

                label dcB = vB->index();

                if (!vB->internalOrBoundaryPoint())
                {
                    dcB = -1;
                }

                label dcOwn = -1;
                label dcNei = -1;

                if (dcA == -1 && dcB == -1)
                {
                    FatalErrorIn("calcDualMesh")
                    << "Attempting to create a face joining "
                        << "two external dual cells "
                        << exit(FatalError);
                }
                else if (dcA == -1 || dcB == -1)
                {
                    // boundary face, find which is the owner

                    if (dcA == -1)
                    {
                        dcOwn = dcB;

                        // reverse face order to correctly orientate normal
                        reverse(newDualFace);
                    }
                    else
                    {
                        dcOwn = dcA;
                    }

                    // Find which patch this face is on by finding the
                    // intersection with the surface of the Delaunay edge
                    // generating the face and identify the region of the
                    // intersection.

                    point ptA = topoint(vA->point());

                    point ptB = topoint(vB->point());

                    pointIndexHit pHit = qSurf_.tree().findLineAny(ptA, ptB);

                    label patchIndex = qSurf_[pHit.index()].region();

                    if (patchIndex == -1)
                    {
                        patchIndex = defaultPatchIndex;

                        WarningIn("Foam::CV3D::calcDualMesh.C")
                            << "Dual face found that is not on a surface "
                            << "patch. Adding to CV3D_default_patch."
                            << endl;
                    }

                    patchFaces[patchIndex].append(newDualFace);
                    patchOwners[patchIndex].append(dcOwn);
                }
                else
                {
                    // internal face, find the lower cell to be the owner

                    if (dcB > dcA)
                    {
                        dcOwn = dcA;
                        dcNei = dcB;
                    }
                    else
                    {
                        dcOwn = dcB;
                        dcNei = dcA;

                        // reverse face order to correctly orientate normal
                        reverse(newDualFace);
                    }

                    faces[dualFacei] = newDualFace;

                    owner[dualFacei] = dcOwn;

                    neighbour[dualFacei] = dcNei;

                    dualFacei++;
                }
            }
            // else
            // {
            //     Info<< verticesOnFace.size()
            //         << " size face not created." << endl;
            // }
        }
    }

    label nInternalFaces = dualFacei;

    faces.setSize(nInternalFaces);

    owner.setSize(nInternalFaces);

    neighbour.setSize(nInternalFaces);

    // ~~~~~~~~ sort owner, reordinging neighbour and faces to match ~~~~~~~~~~~
    // two stage sort for upper triangular order:  sort by owner first, then for
    // each block of owners sort by neighbour

    labelList sortingIndices;

    // Stage 1

    {
        SortableList<label> sortedOwner(owner);

        sortingIndices = sortedOwner.indices();
    }

    {
        labelList copyOwner(owner.size());

        forAll(sortingIndices, sI)
        {
            copyOwner[sI] = owner[sortingIndices[sI]];
        }

        owner = copyOwner;
    }

    {
        labelList copyNeighbour(neighbour.size());

        forAll(sortingIndices, sI)
        {
            copyNeighbour[sI] = neighbour[sortingIndices[sI]];
        }

        neighbour = copyNeighbour;
    }

    {
        faceList copyFaces(faces.size());

        forAll(sortingIndices, sI)
        {
            copyFaces[sI] = faces[sortingIndices[sI]];
        }

        faces = copyFaces;
    }

    // Stage 2

    sortingIndices = -1;

    DynamicList<label> ownerCellJumps;

    // Force first owner entry to be a jump
    ownerCellJumps.append(0);

    for (label o = 1; o < owner.size(); o++)
    {
        if (owner[o] > owner[o-1])
        {
            ownerCellJumps.append(o);
        }
    }

    ownerCellJumps.shrink();

    forAll(ownerCellJumps, oCJ)
    {
        label start = ownerCellJumps[oCJ];

        label length;

        if (oCJ == ownerCellJumps.size() - 1)
        {
            length = owner.size() - start;
        }
        else
        {
            length = ownerCellJumps[oCJ + 1] - start;
        }

        SubList<label> neighbourBlock(neighbour, length, start);

        SortableList<label> sortedNeighbourBlock(neighbourBlock);

        forAll(sortedNeighbourBlock, sNB)
        {
            sortingIndices[start + sNB] =
                sortedNeighbourBlock.indices()[sNB] + start;
        }
    }

    // Perform sort

    {
        labelList copyOwner(owner.size());

        forAll(sortingIndices, sI)
        {
            copyOwner[sI] = owner[sortingIndices[sI]];
        }

        owner = copyOwner;
    }

    {
        labelList copyNeighbour(neighbour.size());

        forAll(sortingIndices, sI)
        {
            copyNeighbour[sI] = neighbour[sortingIndices[sI]];
        }

        neighbour = copyNeighbour;
    }

    {
        faceList copyFaces(faces.size());

        forAll(sortingIndices, sI)
        {
            copyFaces[sI] = faces[sortingIndices[sI]];
        }

        faces = copyFaces;
    }

    // ~~~~~~~~ add patch information ~~~~~~~~~~~

    label nBoundaryFaces = 0;

    forAll(patchFaces, p)
    {
        patchFaces[p].shrink();

        patchOwners[p].shrink();

        patchSizes[p] = patchFaces[p].size();

        patchStarts[p] = nInternalFaces + nBoundaryFaces;

        nBoundaryFaces += patchSizes[p];
    }

    faces.setSize(nInternalFaces + nBoundaryFaces);

    owner.setSize(nInternalFaces + nBoundaryFaces);

    forAll(patchFaces, p)
    {
        forAll(patchFaces[p], f)
        {
            faces[dualFacei] = patchFaces[p][f];

            owner[dualFacei] = patchOwners[p][f];

            dualFacei++;
        }
    }
}


// ************************************************************************* //

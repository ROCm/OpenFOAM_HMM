/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Info << nl << "Calculating Voronoi diagram." << endl;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ dual points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    points.setSize(number_of_cells());

    label dualVerti = 0;

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if
        (
            cit->vertex(0)->internalOrBoundaryPoint()
         || cit->vertex(1)->internalOrBoundaryPoint()
         || cit->vertex(2)->internalOrBoundaryPoint()
         || cit->vertex(3)->internalOrBoundaryPoint()
        )
        {
            cit->cellIndex() = dualVerti;
            points[dualVerti] = topoint(dual(cit));
            dualVerti++;
        }
        else
        {
            cit->cellIndex() = -1;
        }
    }

    points.setSize(dualVerti);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~ dual cell indexing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // resets type and index information for Delaunay vertices
    // assigns an index to the vertices which will be the dual cell index used
    // for owner neighbour assignment

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

    label nPatches = 1;

    patchNames.setSize(nPatches);

    patchNames[0] = "CV3D_default_patch";

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
        if
        (
            eit->first->vertex(eit->second)->internalOrBoundaryPoint()
         || eit->first->vertex(eit->third)->internalOrBoundaryPoint()
        )
        {
            Cell_circulator ccStart = incident_cells(*eit);
            Cell_circulator cc = ccStart;

            DynamicList<label> verticesOnFace;

            do
            {
                if (!is_infinite(cc))
                {
                    if (cc->cellIndex() < 0)
                    {
                        FatalErrorIn("Foam::CV3D::calcDualMesh")
                            << "Dual face uses circumcenter defined by a "
                            << " Delaunay tetrahedron with no internal "
                            << "or boundary points."
                            << exit(FatalError);
                    }

                    verticesOnFace.append(cc->cellIndex());
                }
            } while (++cc != ccStart);

            verticesOnFace.shrink();

            face newDualFace(verticesOnFace);

            Cell_handle c = eit->first;
            Vertex_handle vA = c->vertex(eit->second);
            Vertex_handle vB = c->vertex(eit->third);

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

                // find which patch this face is on;  Hardcoded for now.
                label patchIndex = 0;

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
                    // unsure if CGAL always circulates consistently,
                    // needs to be more rigorous
                    reverse(newDualFace);
                }

                faces[dualFacei] = newDualFace;

                owner[dualFacei] = dcOwn;

                neighbour[dualFacei] = dcNei;

                dualFacei++;
            }
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

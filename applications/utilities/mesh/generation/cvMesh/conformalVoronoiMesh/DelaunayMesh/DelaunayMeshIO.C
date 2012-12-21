/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "DelaunayMesh.H"
#include "fvMesh.H"
#include "pointConversion.H"
#include "wallPolyPatch.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Triangulation>
void Foam::DelaunayMesh<Triangulation>::sortFaces
(
    faceList& faces,
    labelList& owner,
    labelList& neighbour
) const
{
    // Upper triangular order:
    // + owner is sorted in ascending cell order
    // + within each block of equal value for owner, neighbour is sorted in
    //   ascending cell order.
    // + faces sorted to correspond
    // e.g.
    // owner | neighbour
    // 0     | 2
    // 0     | 23
    // 0     | 71
    // 1     | 23
    // 1     | 24
    // 1     | 91

    List<labelPair> ownerNeighbourPair(owner.size());

    forAll(ownerNeighbourPair, oNI)
    {
        ownerNeighbourPair[oNI] = labelPair(owner[oNI], neighbour[oNI]);
    }

    Info<< nl
        << "Sorting faces, owner and neighbour into upper triangular order"
        << endl;

    labelList oldToNew;

    sortedOrder(ownerNeighbourPair, oldToNew);

    oldToNew = invert(oldToNew.size(), oldToNew);

    inplaceReorder(oldToNew, faces);
    inplaceReorder(oldToNew, owner);
    inplaceReorder(oldToNew, neighbour);
}


template<class Triangulation>
void Foam::DelaunayMesh<Triangulation>::addPatches
(
    const label nInternalFaces,
    faceList& faces,
    labelList& owner,
    labelList& patchSizes,
    labelList& patchStarts,
    const List<DynamicList<face> >& patchFaces,
    const List<DynamicList<label> >& patchOwners
) const
{
    label nPatches = patchFaces.size();

    patchSizes.setSize(nPatches, -1);
    patchStarts.setSize(nPatches, -1);

    label nBoundaryFaces = 0;

    forAll(patchFaces, p)
    {
        patchSizes[p] = patchFaces[p].size();
        patchStarts[p] = nInternalFaces + nBoundaryFaces;

        nBoundaryFaces += patchSizes[p];
    }

    faces.setSize(nInternalFaces + nBoundaryFaces);
    owner.setSize(nInternalFaces + nBoundaryFaces);

    label faceI = nInternalFaces;

    forAll(patchFaces, p)
    {
        forAll(patchFaces[p], f)
        {
            faces[faceI] = patchFaces[p][f];
            owner[faceI] = patchOwners[p][f];

            faceI++;
        }
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Triangulation>
void Foam::DelaunayMesh<Triangulation>::printInfo(Ostream& os) const
{
    PrintTable<word, label> triInfoTable("Mesh Statistics");

    triInfoTable.add("Points", Triangulation::number_of_vertices());
    triInfoTable.add("Edges", Triangulation::number_of_finite_edges());
    triInfoTable.add("Faces", Triangulation::number_of_finite_facets());
    triInfoTable.add("Cells", Triangulation::number_of_finite_cells());

    scalar minSize = GREAT;
    scalar maxSize = 0;

    for
    (
        Finite_vertices_iterator vit = Triangulation::finite_vertices_begin();
        vit != Triangulation::finite_vertices_end();
        ++vit
    )
    {
        if (!vit->farPoint())
        {
            minSize = min(vit->targetCellSize(), minSize);
            maxSize = max(vit->targetCellSize(), maxSize);
        }
    }

    Info<< incrIndent;
    triInfoTable.print(Info, true, true);

    Info<< "Size (Min/Max) = "
        << returnReduce(minSize, minOp<scalar>()) << " "
        << returnReduce(maxSize, maxOp<scalar>()) << endl;

    Info<< decrIndent;
}


template<class Triangulation>
Foam::autoPtr<Foam::fvMesh>
Foam::DelaunayMesh<Triangulation>::createMesh
(
    const fileName& name,
    const Time& runTime,
    labelList& vertexMap,
    labelList& cellMap
) const
{
    pointField points(Triangulation::number_of_vertices());
    faceList faces(Triangulation::number_of_finite_facets());
    labelList owner(Triangulation::number_of_finite_facets());
    labelList neighbour(Triangulation::number_of_finite_facets());

    wordList patchNames(1, "cvMesh_defaultPatch");
    wordList patchTypes(1, wallPolyPatch::typeName);

    labelList patchSizes(1, 0);
    labelList patchStarts(1, 0);

    List<DynamicList<face> > patchFaces(1, DynamicList<face>());
    List<DynamicList<label> > patchOwners(1, DynamicList<label>());

    vertexMap.setSize(Triangulation::number_of_vertices());
    cellMap.setSize(Triangulation::number_of_finite_cells());

    // Calculate pts and a map of point index to location in pts.
    label vertI = 0;

    for
    (
        Finite_vertices_iterator vit = Triangulation::finite_vertices_begin();
        vit != Triangulation::finite_vertices_end();
        ++vit
    )
    {
        if (!vit->farPoint())
        {
            vertexMap[vit->index()] = vertI;
            points[vertI] = topoint(vit->point());
            vertI++;
        }
    }

    points.setSize(vertI);

    // Index the cells
    label cellI = 0;

    for
    (
        Finite_cells_iterator cit = Triangulation::finite_cells_begin();
        cit != Triangulation::finite_cells_end();
        ++cit
    )
    {
        if
        (
            !cit->hasFarPoint()
         && !Triangulation::is_infinite(cit)
        )
        {
            cellMap[cit->cellIndex()] = cellI++;
        }
    }

    label faceI = 0;
    labelList verticesOnTriFace(3, -1);
    face newFace(verticesOnTriFace);

    for
    (
        Finite_facets_iterator fit = Triangulation::finite_facets_begin();
        fit != Triangulation::finite_facets_end();
        ++fit
    )
    {
        const Cell_handle c1(fit->first);
        const int oppositeVertex = fit->second;
        const Cell_handle c2(c1->neighbor(oppositeVertex));

        label c1I = Cb::ctFar;
        bool c1Real = false;
        if (!c1->hasFarPoint() && !Triangulation::is_infinite(c1))
        {
            c1I = cellMap[c1->cellIndex()];
            c1Real = true;
        }

        label c2I = Cb::ctFar;
        bool c2Real = false;
        if (!c2->hasFarPoint() && !Triangulation::is_infinite(c2))
        {
            c2I = cellMap[c2->cellIndex()];
            c2Real = true;
        }

        if (!c1Real && !c2Real)
        {
            // Both tets are outside, skip
            continue;
        }

        label ownerCell = -1;
        label neighbourCell = -1;

        for (label i = 0; i < 3; i++)
        {
            verticesOnTriFace[i] = vertexMap
            [
                c1->vertex
                (
                    Triangulation::vertex_triple_index(oppositeVertex, i)
                )->index()
            ];
        }

        newFace = face(verticesOnTriFace);

        if (!c1Real || !c2Real)
        {
            // Boundary face...
            if (!c1Real)
            {
                //... with c1 outside
                ownerCell = c2I;
            }
            else
            {
                // ... with c2 outside
                ownerCell = c1I;

                reverse(newFace);
            }

            patchFaces[0].append(newFace);
            patchOwners[0].append(ownerCell);
        }
        else
        {
            // Internal face...
            if (c1I < c2I)
            {
                // ...with c1 as the ownerCell
                ownerCell = c1I;
                neighbourCell = c2I;

                reverse(newFace);
            }
            else
            {
                // ...with c2 as the ownerCell
                ownerCell = c2I;
                neighbourCell = c1I;
            }

            faces[faceI] = newFace;
            owner[faceI] = ownerCell;
            neighbour[faceI] = neighbourCell;
            faceI++;
        }
    }

    faces.setSize(faceI);
    owner.setSize(faceI);
    neighbour.setSize(faceI);

    sortFaces(faces, owner, neighbour);

    addPatches
    (
        faceI,
        faces,
        owner,
        patchSizes,
        patchStarts,
        patchFaces,
        patchOwners
    );

    autoPtr<fvMesh> meshPtr
    (
        new fvMesh
        (
            IOobject
            (
                name,
                runTime.timeName(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            xferMove(points),
            xferMove(faces),
            xferMove(owner),
            xferMove(neighbour)
        )
    );

    List<polyPatch*> patches(patchStarts.size());

    label nValidPatches = 0;

    forAll(patches, p)
    {
        patches[nValidPatches] = polyPatch::New
        (
            patchTypes[p],
            patchNames[p],
            patchSizes[p],
            patchStarts[p],
            nValidPatches,
            meshPtr().boundaryMesh()
        ).ptr();

        nValidPatches++;
    }

    patches.setSize(nValidPatches);

    meshPtr().addFvPatches(patches);

    return meshPtr;
}


// ************************************************************************* //

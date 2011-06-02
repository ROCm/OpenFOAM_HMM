/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

InClass
    decompositionMethod

\*---------------------------------------------------------------------------*/

#include "decompositionMethod.H"
#include "globalIndex.H"
#include "cyclicPolyPatch.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionMethod, 0);
    defineRunTimeSelectionTable(decompositionMethod, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::decompositionMethod> Foam::decompositionMethod::New
(
    const dictionary& decompositionDict
)
{
    const word methodType(decompositionDict.lookup("method"));

    Info<< "Selecting decompositionMethod " << methodType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(methodType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "decompositionMethod::New"
            "(const dictionary& decompositionDict)"
        )   << "Unknown decompositionMethod "
            << methodType << nl << nl
            << "Valid decompositionMethods are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<decompositionMethod>(cstrIter()(decompositionDict));
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const polyMesh& mesh,
    const pointField& points
)
{
    scalarField weights(points.size(), 1.0);

    return decompose(mesh, points, weights);
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const polyMesh& mesh,
    const labelList& fineToCoarse,
    const pointField& coarsePoints,
    const scalarField& coarseWeights
)
{
    CompactListList<label> coarseCellCells;
    calcCellCells
    (
        mesh,
        fineToCoarse,
        coarsePoints.size(),
        coarseCellCells
    );

    // Decompose based on agglomerated points
    labelList coarseDistribution
    (
        decompose
        (
            coarseCellCells(),
            coarsePoints,
            coarseWeights
        )
    );

    // Rework back into decomposition for original mesh_
    labelList fineDistribution(fineToCoarse.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = coarseDistribution[fineToCoarse[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const polyMesh& mesh,
    const labelList& fineToCoarse,
    const pointField& coarsePoints
)
{
    scalarField cWeights(coarsePoints.size(), 1.0);

    return decompose
    (
        mesh,
        fineToCoarse,
        coarsePoints,
        cWeights
    );
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const labelListList& globalCellCells,
    const pointField& cc
)
{
    scalarField cWeights(cc.size(), 1.0);

    return decompose(globalCellCells, cc, cWeights);
}


// Return the minimum face between two cells. Only relevant for
// cells with multiple faces inbetween.
Foam::label Foam::decompositionMethod::masterFace
(
    const polyMesh& mesh,
    const label own,
    const label nei
)
{
    label minFaceI = labelMax;

    // Count multiple faces between own and nei only once
    const cell& ownFaces = mesh.cells()[own];
    forAll(ownFaces, i)
    {
        label otherFaceI = ownFaces[i];

        if (mesh.isInternalFace(otherFaceI))
        {
            label nbrCellI =
            (
                mesh.faceNeighbour()[otherFaceI] != own
              ? mesh.faceNeighbour()[otherFaceI]
              : mesh.faceOwner()[otherFaceI]
            );

            if (nbrCellI == nei)
            {
                minFaceI = min(minFaceI, otherFaceI);
            }
        }
    }

    return minFaceI;
}


//void Foam::decompositionMethod::calcCSR
//(
//    const polyMesh& mesh,
//    List<int>& adjncy,
//    List<int>& xadj
//)
//{
//    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
//
//    // Make Metis CSR (Compressed Storage Format) storage
//    //   adjncy      : contains neighbours (= edges in graph)
//    //   xadj(celli) : start of information in adjncy for celli
//
//
//    // Count unique faces between cells
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//    labelList nFacesPerCell(mesh.nCells(), 0);
//
//    // Internal faces
//    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
//    {
//        label own = mesh.faceOwner()[faceI];
//        label nei = mesh.faceNeighbour()[faceI];
//
//        if (faceI == masterFace(mesh, own, nei))
//        {
//            nFacesPerCell[own]++;
//            nFacesPerCell[nei]++;
//        }
//    }
//
//    // Coupled faces. Only cyclics done.
//    HashSet<edge, Hash<edge> > cellPair(mesh.nFaces()-mesh.nInternalFaces());
//
//    forAll(pbm, patchI)
//    {
//        if
//        (
//            isA<cyclicPolyPatch>(pbm[patchI])
//         && refCast<const cyclicPolyPatch>(pbm[patchI]).owner()
//        )
//        {
//            const cyclicPolyPatch& cycPatch = refCast<const cyclicPolyPatch>
//            (
//                pbm[patchI]
//            );
//
//            const labelUList& faceCells = cycPatch.faceCells();
//            const labelUList& nbrCells =
//                cycPatch.neighbPatch().faceCells();
//
//            forAll(faceCells, facei)
//            {
//                label own = faceCells[facei];
//                label nei = nbrCells[facei];
//
//                if (cellPair.insert(edge(own, nei)))
//                {
//                    nFacesPerCell[own]++;
//                    nFacesPerCell[nei]++;
//                }
//            }
//        }
//    }
//
//
//    // Size tables
//    // ~~~~~~~~~~~
//
//    // Sum nFacesPerCell
//    xadj.setSize(mesh.nCells()+1);
//
//    label nConnections = 0;
//
//    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
//    {
//        xadj[cellI] = nConnections;
//        nConnections += nFacesPerCell[cellI];
//    }
//    xadj[mesh.nCells()] = nConnections;
//    adjncy.setSize(nConnections);
//
//
//
//    // Fill tables
//    // ~~~~~~~~~~~
//
//    nFacesPerCell = 0;
//
//    // Internal faces
//    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
//    {
//        label own = mesh.faceOwner()[faceI];
//        label nei = mesh.faceNeighbour()[faceI];
//
//        if (faceI == masterFace(mesh, own, nei))
//        {
//            adjncy[xadj[own] + nFacesPerCell[own]++] = nei;
//            adjncy[xadj[nei] + nFacesPerCell[nei]++] = own;
//        }
//    }
//
//    // Coupled faces. Only cyclics done.
//    cellPair.clear();
//    forAll(pbm, patchI)
//    {
//        if
//        (
//            isA<cyclicPolyPatch>(pbm[patchI])
//         && refCast<const cyclicPolyPatch>(pbm[patchI]).owner()
//        )
//        {
//            const cyclicPolyPatch& cycPatch = refCast<const cyclicPolyPatch>
//            (
//                pbm[patchI]
//            );
//
//            const labelUList& faceCells = cycPatch.faceCells();
//            const labelUList& nbrCells =
//                cycPatch.neighbPatch().faceCells();
//
//            forAll(faceCells, facei)
//            {
//                label own = faceCells[facei];
//                label nei = nbrCells[facei];
//
//                if (cellPair.insert(edge(own, nei)))
//                {
//                    adjncy[xadj[own] + nFacesPerCell[own]++] = nei;
//                    adjncy[xadj[nei] + nFacesPerCell[nei]++] = own;
//                }
//            }
//        }
//    }
//}
//
//
//// From cell-cell connections to Metis format (like CompactListList)
//void Foam::decompositionMethod::calcCSR
//(
//    const labelListList& cellCells,
//    List<int>& adjncy,
//    List<int>& xadj
//)
//{
//    labelHashSet nbrCells;
//
//    // Count number of internal faces
//    label nConnections = 0;
//
//    forAll(cellCells, coarseI)
//    {
//        nbrCells.clear();
//
//        const labelList& cCells = cellCells[coarseI];
//
//        forAll(cCells, i)
//        {
//            if (nbrCells.insert(cCells[i]))
//            {
//                nConnections++;
//            }
//        }
//    }
//
//    // Create the adjncy array as twice the size of the total number of
//    // internal faces
//    adjncy.setSize(nConnections);
//
//    xadj.setSize(cellCells.size()+1);
//
//
//    // Fill in xadj
//    // ~~~~~~~~~~~~
//    label freeAdj = 0;
//
//    forAll(cellCells, coarseI)
//    {
//        xadj[coarseI] = freeAdj;
//
//        nbrCells.clear();
//
//        const labelList& cCells = cellCells[coarseI];
//
//        forAll(cCells, i)
//        {
//            if (nbrCells.insert(cCells[i]))
//            {
//                adjncy[freeAdj++] = cCells[i];
//            }
//        }
//    }
//    xadj[cellCells.size()] = freeAdj;
//}


void Foam::decompositionMethod::calcCellCells
(
    const polyMesh& mesh,
    const labelList& agglom,
    const label nCoarse,
    CompactListList<label>& cellCells
)
{
    // Create global cell numbers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nAgglom;
    if (agglom.size())
    {
        nAgglom = max(agglom)+1;
    }
    else
    {
        nAgglom = 0;
    }
    globalIndex globalAgglom(nAgglom);

    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Get agglomerate owner on other side of coupled faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList globalNeighbour(mesh.nFaces()-mesh.nInternalFaces());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            label bFaceI = pp.start() - mesh.nInternalFaces();

            forAll(pp, i)
            {
                globalNeighbour[bFaceI] = globalAgglom.toGlobal
                (
                    agglom[faceOwner[faceI]]
                );

                bFaceI++;
                faceI++;
            }
        }
    }

    // Get the cell on the other side of coupled patches
    syncTools::swapBoundaryFaceList(mesh, globalNeighbour);


    // Count number of faces (internal + coupled)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Number of faces per coarse cell
    labelList nFacesPerCell(nCoarse, 0);

    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        label own = agglom[faceOwner[faceI]];
        label nei = agglom[faceNeighbour[faceI]];

        nFacesPerCell[own]++;
        nFacesPerCell[nei]++;
    }

    // Handle coupled faces. In case of agglomeration you might end up
    // with multiple connections ('faces') between the same two agglomerations.
    // This is illegal so mark.
    HashSet<labelPair, Hash<labelPair> > cellPair
    (
        mesh.nFaces()-mesh.nInternalFaces()
    );

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            label bFaceI = pp.start()-mesh.nInternalFaces();

            forAll(pp, i)
            {
                label own = agglom[faceOwner[faceI]];
                label globalNei = globalNeighbour[bFaceI];
                if
                (
                    globalAgglom.toGlobal(own) != globalNei
                 && cellPair.insert(labelPair(own, globalNei))
                )
                {
                    nFacesPerCell[own]++;
                }
                faceI++;
                bFaceI++;
            }
        }
    }


    // Fill in offset and data
    // ~~~~~~~~~~~~~~~~~~~~~~~

    cellCells.setSize(nFacesPerCell);

    nFacesPerCell = 0;

    labelList& m = cellCells.m();
    const labelList& offsets = cellCells.offsets();

    // For internal faces is just offsetted owner and neighbour
    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        label own = agglom[faceOwner[faceI]];
        label nei = agglom[faceNeighbour[faceI]];

        m[offsets[own] + nFacesPerCell[own]++] = globalAgglom.toGlobal(nei);
        m[offsets[nei] + nFacesPerCell[nei]++] = globalAgglom.toGlobal(own);
    }

    // For boundary faces is offsetted coupled neighbour
    cellPair.clear();
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            label bFaceI = pp.start()-mesh.nInternalFaces();

            forAll(pp, i)
            {
                label own = agglom[faceOwner[faceI]];
                label globalNei = globalNeighbour[bFaceI];
                if (cellPair.insert(labelPair(own, globalNei)))
                {
                    m[offsets[own] + nFacesPerCell[own]++] = globalNei;
                }
                faceI++;
                bFaceI++;
            }
        }
    }


    // Check for duplicates connections between cells as a postprocessing step
    // (since quite rare)
    label startIndex = 0;
    label newIndex = 0;
    labelHashSet nbrCells;
    forAll(cellCells, cellI)
    {
        nbrCells.clear();

        label& endIndex = cellCells.offsets()[cellI+1];

        for (label i = startIndex; i < endIndex; i++)
        {
            if (nbrCells.insert(cellCells.m()[i]))
            {
                cellCells.m()[newIndex++] = cellCells.m()[i];
            }
        }

        startIndex = endIndex;
        endIndex = newIndex;
    }
}


// ************************************************************************* //

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

\*---------------------------------------------------------------------------*/

#include "extendedStencil.H"
#include "globalIndex.H"
#include "syncTools.H"
#include "SortableList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculates per face a list of global cell/face indices.
void Foam::extendedStencil::calcFaceStencils
(
    const polyMesh& mesh,
    const globalIndex& globalNumbering
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const label nBnd = mesh.nFaces()-mesh.nInternalFaces();
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();


    // Determine neighbouring global cell or boundary face
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList neiGlobal(nBnd);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        label faceI = pp.start();

        if (pp.coupled())
        {
            // For coupled faces get the cell on the other side
            forAll(pp, i)
            {
                label bFaceI = faceI-mesh.nInternalFaces();
                neiGlobal[bFaceI] = globalNumbering.toGlobal(own[faceI]);
                faceI++;
            }
        }
        else if (isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label bFaceI = faceI-mesh.nInternalFaces();
                neiGlobal[bFaceI] = -1;
                faceI++;
            }
        }
        else
        {
            // For noncoupled faces get the boundary face.
            forAll(pp, i)
            {
                label bFaceI = faceI-mesh.nInternalFaces();
                neiGlobal[bFaceI] =
                    globalNumbering.toGlobal(mesh.nCells()+bFaceI);
                faceI++;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh, neiGlobal, false);


    // Determine cellCells in global numbering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList globalCellCells(mesh.nCells());
    forAll(globalCellCells, cellI)
    {
        const cell& cFaces = mesh.cells()[cellI];

        labelList& cCells = globalCellCells[cellI];

        cCells.setSize(cFaces.size());

        // Collect neighbouring cells/faces
        label nNbr = 0;
        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            if (mesh.isInternalFace(faceI))
            {
                label nbrCellI = own[faceI];
                if (nbrCellI == cellI)
                {
                    nbrCellI = nei[faceI];
                }
                cCells[nNbr++] = globalNumbering.toGlobal(nbrCellI);
            }
            else
            {
                label nbrCellI = neiGlobal[faceI-mesh.nInternalFaces()];
                if (nbrCellI != -1)
                {
                    cCells[nNbr++] = nbrCellI;
                }
            }
        }
        cCells.setSize(nNbr);
    }


    // Determine neighbouring global cell Cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList neiGlobalCellCells(nBnd);
    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        neiGlobalCellCells[faceI-mesh.nInternalFaces()] =
            globalCellCells[own[faceI]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiGlobalCellCells, false);



    // Construct stencil in global numbering
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    stencil_.setSize(mesh.nFaces());

    labelHashSet faceStencil;

    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        faceStencil.clear();
        label globalOwn = globalNumbering.toGlobal(own[faceI]);
        faceStencil.insert(globalOwn);
        const labelList& ownCCells = globalCellCells[own[faceI]];
        forAll(ownCCells, i)
        {
            faceStencil.insert(ownCCells[i]);
        }

        label globalNei = globalNumbering.toGlobal(nei[faceI]);
        faceStencil.insert(globalNei);
        const labelList& neiCCells = globalCellCells[nei[faceI]];
        forAll(neiCCells, i)
        {
            faceStencil.insert(neiCCells[i]);
        }

        // Guarantee owner first, neighbour second.
        stencil_[faceI].setSize(faceStencil.size());
        label n = 0;
        stencil_[faceI][n++] = globalOwn;
        stencil_[faceI][n++] = globalNei;
        forAllConstIter(labelHashSet, faceStencil, iter)
        {
            if (iter.key() != globalOwn && iter.key() != globalNei)
            {
                stencil_[faceI][n++] = iter.key();
            }
        }
        //Pout<< "internalface:" << faceI << " toc:" << faceStencil.toc()
        //    << " stencil:" << stencil_[faceI] << endl;
    }
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        label faceI = pp.start();

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                faceStencil.clear();
                label globalOwn = globalNumbering.toGlobal(own[faceI]);
                faceStencil.insert(globalOwn);
                const labelList& ownCCells = globalCellCells[own[faceI]];
                forAll(ownCCells, i)
                {
                    faceStencil.insert(ownCCells[i]);
                }
                // Get the coupled cell
                label globalNei = neiGlobal[faceI-mesh.nInternalFaces()];
                faceStencil.insert(globalNei);
                // And the neighbours of the coupled cell
                const labelList& neiCCells =
                    neiGlobalCellCells[faceI-mesh.nInternalFaces()];
                forAll(neiCCells, i)
                {
                    faceStencil.insert(neiCCells[i]);
                }

                // Guarantee owner first, neighbour second.
                stencil_[faceI].setSize(faceStencil.size());
                label n = 0;
                stencil_[faceI][n++] = globalOwn;
                stencil_[faceI][n++] = globalNei;
                forAllConstIter(labelHashSet, faceStencil, iter)
                {
                    if (iter.key() != globalOwn && iter.key() != globalNei)
                    {
                        stencil_[faceI][n++] = iter.key();
                    }
                }

                //Pout<< "coupledface:" << faceI
                //    << " toc:" << faceStencil.toc()
                //    << " stencil:" << stencil_[faceI] << endl;

                faceI++;
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                faceStencil.clear();
                label globalOwn = globalNumbering.toGlobal(own[faceI]);
                faceStencil.insert(globalOwn);
                const labelList& ownCCells = globalCellCells[own[faceI]];
                forAll(ownCCells, i)
                {
                    faceStencil.insert(ownCCells[i]);
                }


                // Guarantee owner first, neighbour second.
                stencil_[faceI].setSize(faceStencil.size());
                label n = 0;
                stencil_[faceI][n++] = globalOwn;
                forAllConstIter(labelHashSet, faceStencil, iter)
                {
                    if (iter.key() != globalOwn)
                    {
                        stencil_[faceI][n++] = iter.key();
                    }
                }

                //Pout<< "boundaryface:" << faceI
                //    << " toc:" << faceStencil.toc()
                //    << " stencil:" << stencil_[faceI] << endl;

                faceI++;
            }
        }
    }
}


// Calculates extended stencil. This is per face
// - owner
// - cellCells of owner
// - neighbour
// - cellCells of neighbour
// It comes in two parts:
// - a map which collects/distributes all necessary data in a compact array
// - the stencil (a labelList per face) which is a set of indices into this
//   compact array.
// The compact array is laid out as follows:
// - first data for current processor (Pstream::myProcNo())
//      - all cells
//      - all boundary faces
// - then per processor
//      - all used cells and boundary faces
void Foam::extendedStencil::calcExtendedFaceStencil(const polyMesh& mesh)
{
    const label nBnd = mesh.nFaces()-mesh.nInternalFaces();

    // Global numbering for cells and boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    globalIndex globalNumbering(mesh.nCells()+nBnd);


    // Calculate stencil in global cell indices
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    calcFaceStencils(mesh, globalNumbering);


    // Convert stencil to schedule
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // We now know what information we need from other processors. This needs
    // to be converted into what information I need to send as well
    // (mapDistribute)


    // 1. Construct per processor compact addressing of the global cells
    //    needed. The ones from the local processor are not included since
    //    these are always all needed.
    List<Map<label> > globalToProc(Pstream::nProcs());
    {
        const labelList& procPatchMap = mesh.globalData().procPatchMap();
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Presize with (as estimate) size of patch to neighbour.
        forAll(procPatchMap, procI)
        {
            if (procPatchMap[procI] != -1)
            {
                globalToProc[procI].resize
                (
                    patches[procPatchMap[procI]].size()
                );
            }
        }

        // Collect all (non-local) globalcells/faces needed.
        forAll(stencil_, faceI)
        {
            const labelList& stencilCells = stencil_[faceI];

            forAll(stencilCells, i)
            {
                label globalCellI = stencilCells[i];
                label procI = globalNumbering.whichProcID(stencilCells[i]);

                if (procI != Pstream::myProcNo())
                {
                    label nCompact = globalToProc[procI].size();
                    globalToProc[procI].insert(globalCellI, nCompact);
                }
            }
        }
        // Sort global cells needed (not really necessary)
        forAll(globalToProc, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                Map<label>& globalMap = globalToProc[procI];

                SortableList<label> sorted(globalMap.toc());

                forAll(sorted, i)
                {
                    Map<label>::iterator iter = globalMap.find(sorted[i]);
                    iter() = i;
                }
            }
        }


        // forAll(globalToProc, procI)
        // {
        //     Pout<< "From processor:" << procI << " want cells/faces:" << endl;
        //     forAllConstIter(Map<label>, globalToProc[procI], iter)
        //     {
        //         Pout<< "    global:" << iter.key()
        //             << " local:" << globalNumbering.toLocal(procI, iter.key())
        //             << endl;
        //     }
        //     Pout<< endl;
        // }
    }


    // 2. The overall compact addressing is
    // - myProcNo first
    // - all other processors consecutively

    labelList compactStart(Pstream::nProcs());
    compactStart[Pstream::myProcNo()] = 0;
    label nCompact = mesh.nCells()+nBnd;
    forAll(compactStart, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            compactStart[procI] = nCompact;
            nCompact += globalToProc[procI].size();

            // Pout<< "Data wanted from " << procI << " starts at "
            //     << compactStart[procI] << endl;
        }
    }
    // Pout<< "Overall cells needed:" << nCompact << endl;


    // 3. Find out what to receive/send in compact addressing.
    labelListList recvCompact(Pstream::nProcs());
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            labelList wantedGlobals(globalToProc[procI].size());
            recvCompact[procI].setSize(globalToProc[procI].size());

            label i = 0;
            forAllConstIter(Map<label>, globalToProc[procI], iter)
            {
                wantedGlobals[i] = iter.key();
                recvCompact[procI][i] = compactStart[procI]+iter();
                i++;
            }

            // Pout<< "From proc:" << procI
            //     << " I need (globalcells):" << wantedGlobals
            //     << " which are my compact:" << recvCompact[procI]
            //     << endl;

            // Send the global cell numbers I need from procI
            OPstream str(Pstream::blocking, procI);
            str << wantedGlobals;
        }
        else
        {
            recvCompact[procI] =
                compactStart[procI]
              + identity(mesh.nCells()+nBnd);
        }
    }
    labelListList sendCompact(Pstream::nProcs());
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            // See what neighbour wants to receive (= what I need to send)

            IPstream str(Pstream::blocking, procI);
            labelList globalCells(str);

            labelList& procCompact = sendCompact[procI];
            procCompact.setSize(globalCells.size());

            // Convert from globalCells (all on my processor!) into compact
            // addressing
            forAll(globalCells, i)
            {
                label cellI = globalNumbering.toLocal(globalCells[i]);
                procCompact[i] = compactStart[Pstream::myProcNo()]+cellI;
            }
        }
        else
        {
            sendCompact[procI] = recvCompact[procI];
        }
    }

    // Convert stencil to compact numbering
    forAll(stencil_, faceI)
    {
        labelList& stencilCells = stencil_[faceI];

        forAll(stencilCells, i)
        {
            label globalCellI = stencilCells[i];
            label procI = globalNumbering.whichProcID(globalCellI);
            if (procI != Pstream::myProcNo())
            {
                label localCompact = globalToProc[procI][globalCellI];
                stencilCells[i] = compactStart[procI]+localCompact;
            }
            else
            {
                label localCompact = globalNumbering.toLocal(globalCellI);
                stencilCells[i] = compactStart[procI]+localCompact;
            }

        }
    }
    // Pout<< "***stencil_:" << stencil_ << endl;

    // Constuct map for distribution of compact data.
    mapPtr_.reset
    (
        new mapDistribute
        (
            nCompact,
            sendCompact,
            recvCompact,
            true            // reuse send/recv maps.
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedStencil::extendedStencil
(
    const mapDistribute& map,
    const labelListList& stencil
)
:
    mapPtr_
    (
        autoPtr<mapDistribute>
        (
            new mapDistribute
            (
                map.constructSize(),
                map.subMap(),
                map.constructMap()
            )
        )
    ),
    stencil_(stencil)
{}


Foam::extendedStencil::extendedStencil(const polyMesh& mesh)
{
    calcExtendedFaceStencil(mesh);
}


// ************************************************************************* //

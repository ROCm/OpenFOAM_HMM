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

#include "syncTools.H"
#include "parMetisDecomp.H"
#include "metisDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "polyMesh.H"
#include "Time.H"
#include "labelIOField.H"

#include <mpi.h>

extern "C"
{
#   include "parmetis.h"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(parMetisDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        parMetisDecomp,
        dictionaryMesh
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parMetisDecomp::parMetisDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::parMetisDecomp::decompose(const pointField& points)
{
    // For running sequential ...
    if (Pstream::nProcs() <= 1)
    {
        return metisDecomp(decompositionDict_, mesh_).decompose(points);
    }

    //
    // Make Metis Distributed CSR (Compressed Storage Format) storage
    //   adjncy      : contains cellCells (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    //

    // Create global cell numbers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Get number of cells on all processors
    labelList nLocalCells(Pstream::nProcs());
    nLocalCells[Pstream::myProcNo()] = mesh_.nCells();
    Pstream::gatherList(nLocalCells);
    Pstream::scatterList(nLocalCells);

    // Get cell offsets.
    labelList cellOffsets(Pstream::nProcs()+1);
    label nGlobalCells = 0;
    forAll(nLocalCells, procI)
    {
        cellOffsets[procI] = nGlobalCells;
        nGlobalCells += nLocalCells[procI];
    }
    cellOffsets[Pstream::nProcs()] = nGlobalCells;

    label myOffset = cellOffsets[Pstream::myProcNo()];



    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();


    // Get renumbered owner on other side of coupled faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList globalNeighbour(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            label bFaceI = pp.start() - mesh_.nInternalFaces();

            forAll(pp, i)
            {
                globalNeighbour[bFaceI++] = faceOwner[faceI++] + myOffset;
            }
        }
    }

    // Get the cell on the other side of coupled patches
    syncTools::swapBoundaryFaceList(mesh_, globalNeighbour, false);


    // Count number of faces (internal + coupled)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Number of faces per cell
    labelList nFacesPerCell(mesh_.nCells(), 0);

    // Number of coupled faces
    label nCoupledFaces = 0;

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        nFacesPerCell[faceOwner[faceI]]++;
        nFacesPerCell[faceNeighbour[faceI]]++;
    }
    // Handle coupled faces
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                nCoupledFaces++;
                nFacesPerCell[faceOwner[faceI++]]++;
            }
        }
    }            


    // Fill in xadj
    // ~~~~~~~~~~~~

    labelField xadj(mesh_.nCells()+1, -1);

    label freeAdj = 0;

    for (label cellI = 0; cellI < mesh_.nCells(); cellI++)
    {
        xadj[cellI] = freeAdj;

        freeAdj += nFacesPerCell[cellI];
    }
    xadj[mesh_.nCells()] = freeAdj;



    // Fill in adjncy
    // ~~~~~~~~~~~~~~

    labelField adjncy(2*mesh_.nInternalFaces() + nCoupledFaces, -1);

    nFacesPerCell = 0;

    // For internal faces is just offsetted owner and neighbour
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = faceOwner[faceI];
        label nei = faceNeighbour[faceI];

        adjncy[xadj[own] + nFacesPerCell[own]++] = nei + myOffset;
        adjncy[xadj[nei] + nFacesPerCell[nei]++] = own + myOffset;
    }
    // For boundary faces is offsetted coupled neighbour
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            label bFaceI = pp.start()-mesh_.nInternalFaces();

            forAll(pp, i)
            {
                label own = faceOwner[faceI];
                adjncy[xadj[own] + nFacesPerCell[own]++] =
                    globalNeighbour[bFaceI];

                faceI++;
                bFaceI++;
            }
        }
    }            



    // C style numbering
    int numFlag = 0;

    // Number of dimensions
    int nDims = 3;

    // cell centres
    Field<floatScalar> xyz(nDims*mesh_.nCells());
    const pointField& cellCentres = mesh_.cellCentres();
    label compI = 0;
    forAll(cellCentres, cellI)
    {
        const point& cc = cellCentres[cellI];
        xyz[compI++] = float(cc.x());
        xyz[compI++] = float(cc.y());
        xyz[compI++] = float(cc.z());
    }


    // decomposition options. 0 = use defaults
    labelList options(3, 0);
    //options[0] = 1;     // don't use defaults but use values below
    //options[1] = -1;    // full debug info
    //options[2] = 15;    // random number seed


    // cell weights (so on the vertices of the dual)
    labelField cellWeights;

    // face weights (so on the edges of the dual)
    labelField faceWeights;

    // Check for user supplied weights and decomp options
    if (decompositionDict_.found("metisCoeffs"))
    {
        dictionary parMetisDecompCoeffs
        (
            decompositionDict_.subDict("metisCoeffs")
        );

        if (parMetisDecompCoeffs.found("cellWeightsFile"))
        {
            word cellWeightsFile
            (
                parMetisDecompCoeffs.lookup("cellWeightsFile")
            );

            Info<< "parMetisDecomp : Using cell-based weights read from "
                << cellWeightsFile << endl;

            labelIOField cellIOWeights
            (
                IOobject
                (
                    cellWeightsFile,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
            cellWeights.transfer(cellIOWeights);

            if (cellWeights.size() != mesh_.nCells())
            {
                FatalErrorIn("parMetisDecomp::decompose(const pointField&)")
                    << "Number of cell weights " << cellWeights.size()
                    << " read from " << cellIOWeights.objectPath()
                    << " does not equal number of cells " << mesh_.nCells()
                    << exit(FatalError);
            }
        }

        if (parMetisDecompCoeffs.found("faceWeightsFile"))
        {
            word faceWeightsFile
            (
                parMetisDecompCoeffs.lookup("faceWeightsFile")
            );

            Info<< "parMetisDecomp : Using face-based weights read from "
                << faceWeightsFile << endl;

            labelIOField weights
            (
                IOobject
                (
                    faceWeightsFile,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );

            if (weights.size() != mesh_.nFaces())
            {
                FatalErrorIn("parMetisDecomp::decompose(const pointField&)")
                    << "Number of face weights " << weights.size()
                    << " does not equal number of internal and boundary faces "
                    << mesh_.nFaces()
                    << exit(FatalError);
            }

            faceWeights.setSize(2*mesh_.nInternalFaces()+nCoupledFaces);

            // Assume symmetric weights. Keep same ordering as adjncy.
            nFacesPerCell = 0;

            // Handle internal faces
            for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
            {
                label w = weights[faceI];

                label own = faceOwner[faceI];
                label nei = faceNeighbour[faceI];

                faceWeights[xadj[own] + nFacesPerCell[own]++] = w;
                faceWeights[xadj[nei] + nFacesPerCell[nei]++] = w;
            }
            // Coupled boundary faces
            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];

                if (pp.coupled())
                {
                    label faceI = pp.start();

                    forAll(pp, i)
                    {
                        label w = weights[faceI];
                        label own = faceOwner[faceI];
                        adjncy[xadj[own] + nFacesPerCell[own]++] = w;
                        faceI++;
                    }
                }
            }            
        }

        if (parMetisDecompCoeffs.found("options"))
        {
            parMetisDecompCoeffs.lookup("options") >> options;

            Info<< "Using Metis options     " << options
                << endl << endl;

            if (options.size() != 3)
            {
                FatalErrorIn("parMetisDecomp::decompose(const pointField&)")
                    << "Number of options " << options.size()
                    << " should be three." << exit(FatalError);
            }
        }
    }



    // Make sure every domain has at least one cell
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (Metis falls over with zero sized domains)
    // Trickle cells from processors that have them down to those that
    // don't.


    // Number of cells to send down (is same as number of cells next processor
    // has to receive)
    labelList nSendCells(Pstream::nProcs(), 0);

    for (label procI = nLocalCells.size()-1; procI >=1; procI--)
    {
        if (nLocalCells[procI]-nSendCells[procI] < 1)
        {
            nSendCells[procI-1] = nSendCells[procI]-nLocalCells[procI]+1;
        }
    }

    // First receive (so increasing the sizes of all arrays)

    if (Pstream::myProcNo() >= 1 && nSendCells[Pstream::myProcNo()-1] > 0)
    {
        // Receive cells from previous processor
        IPstream fromPrevProc(Pstream::blocking, Pstream::myProcNo()-1);

        labelField prevXadj(fromPrevProc);
        labelField prevAdjncy(fromPrevProc);
        Field<floatScalar> prevXyz(fromPrevProc);
        labelField prevCellWeights(fromPrevProc);
        labelField prevFaceWeights(fromPrevProc);

        // Insert adjncy
        prepend(prevAdjncy, adjncy);
        // Adapt offsets and prepend xadj
        xadj += prevAdjncy.size();
        prepend(prevXadj, xadj);
        // Coords
        prepend(prevXyz, xyz);
        // Weights
        prepend(prevCellWeights, cellWeights);
        prepend(prevFaceWeights, faceWeights);
    }


    // Send to my next processor

    if (nSendCells[Pstream::myProcNo()] > 0)
    {
        // Send cells to next processor
        OPstream toNextProc(Pstream::blocking, Pstream::myProcNo()+1);

        label nCells = nSendCells[Pstream::myProcNo()];
        label startCell = xadj.size()-1 - nCells;
        label startFace = xadj[startCell];
        label nFaces = adjncy.size()-startFace;

        // Send for all cell data: last nCells elements
        // Send for all face data: last nFaces elements
        toNextProc
            << labelField::subField(xadj, nCells, startCell)-startFace
            << labelField::subField(adjncy, nFaces, startFace)
            << SubField<floatScalar>(xyz, nDims*nCells, nDims*startCell)
            <<
            (
                (cellWeights.size() > 0)
              ? static_cast<const labelField&>
                (
                    labelField::subField(cellWeights, nCells, startCell)
                )
              : labelField(0)
            )
            <<
            (
                (faceWeights.size() > 0)
              ? static_cast<const labelField&>
                (
                    labelField::subField(faceWeights, nFaces, startFace)
                )
              : labelField(0)
            );

        // Remove data that has been sent
        if (faceWeights.size() > 0)
        {
            faceWeights.setSize(faceWeights.size()-nFaces);
        }
        if (cellWeights.size() > 0)
        {
            cellWeights.setSize(cellWeights.size()-nCells);
        }
        xyz.setSize(xyz.size()-nDims*nCells);
        adjncy.setSize(adjncy.size()-nFaces);
        xadj.setSize(xadj.size() - nCells);
    }



    // Adapt number of cells
    forAll(nSendCells, procI)
    {
        // Sent cells
        nLocalCells[procI] -= nSendCells[procI];

        if (procI >= 1)
        {
            // Received cells
            nLocalCells[procI] += nSendCells[procI-1];
        }
    }
    // Adapt cellOffsets
    nGlobalCells = 0;
    forAll(nLocalCells, procI)
    {
        cellOffsets[procI] = nGlobalCells;
        nGlobalCells += nLocalCells[procI];
    }


    // Weight info
    int wgtFlag = 0;
    label* vwgtPtr = NULL;
    label* adjwgtPtr = NULL;

    if (cellWeights.size() > 0)
    {
        vwgtPtr = cellWeights.begin();
        wgtFlag += 2;       // Weights on vertices
    }
    if (faceWeights.size() > 0)
    {
        adjwgtPtr = faceWeights.begin();
        wgtFlag += 1;       // Weights on edges
    }


    // Number of weights or balance constraints
    int nCon = 1;
    // Per processor, per constraint the weight
    Field<floatScalar> tpwgts(nCon*nProcessors_, 1./nProcessors_);
    // Imbalance tolerance
    Field<floatScalar> ubvec(nCon, 1.02);
    if (nProcessors_ == 1)
    {
        // If only one processor there is no imbalance.
        ubvec[0] = 1;
    }

    MPI_Comm comm = MPI_COMM_WORLD;

    // output: cell -> processor addressing
    labelList finalDecomp(nLocalCells[Pstream::myProcNo()]);

    // output: number of cut edges
    int edgeCut = 0;


    ParMETIS_V3_PartGeomKway
    (
        cellOffsets.begin(),    // vtxDist
        xadj.begin(),
        adjncy.begin(),
        vwgtPtr,                // vertexweights
        adjwgtPtr,              // edgeweights
        &wgtFlag,
        &numFlag,
        &nDims,
        xyz.begin(),
        &nCon,
        &nProcessors_,          // nParts
        tpwgts.begin(),
        ubvec.begin(),
        options.begin(),
        &edgeCut,
        finalDecomp.begin(),
        &comm
    );


    // If we sent cells across make sure we undo it
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Receive back from next processor if I sent something
    if (nSendCells[Pstream::myProcNo()] > 0)
    {
        IPstream fromNextProc(Pstream::blocking, Pstream::myProcNo()+1);

        labelList nextFinalDecomp(fromNextProc);

        append(nextFinalDecomp, finalDecomp);
    }

    // Send back to previous processor.
    if (Pstream::myProcNo() >= 1 && nSendCells[Pstream::myProcNo()-1] > 0)
    {
        OPstream toPrevProc(Pstream::blocking, Pstream::myProcNo()-1);

        label nToPrevious = nSendCells[Pstream::myProcNo()-1];

        toPrevProc <<
            SubList<label>
            (
                finalDecomp,
                nToPrevious,
                finalDecomp.size()-nToPrevious
            );

        // Remove locally what has been sent
        finalDecomp.setSize(finalDecomp.size()-nToPrevious);
    }

    return finalDecomp;
}


// ************************************************************************* //

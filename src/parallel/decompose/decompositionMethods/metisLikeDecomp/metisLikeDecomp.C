/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "metisLikeDecomp.H"
#include "Time.H"
#include "globalIndex.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::metisLikeDecomp::decomposeGeneral
(
    const labelList& adjncy,
    const labelList& xadj,
    const List<scalar>& cWeights,
    labelList& decomp
) const
{
    if (!Pstream::parRun())
    {
        return decomposeSerial
        (
            adjncy,
            xadj,
            cWeights,
            decomp
        );
    }

    if (debug)
    {
        Info<< type() << "Decomp : running in parallel."
            << " Decomposing all of graph on master processor." << endl;
    }
    globalIndex globalCells(xadj.size()-1);
    label nTotalConnections = returnReduce(adjncy.size(), sumOp<label>());

    // Send all to master. Use scheduled to save some storage.
    if (Pstream::master())
    {
        List<label> allAdjncy(nTotalConnections);
        List<label> allXadj(globalCells.size()+1);
        List<scalar> allWeights(globalCells.size());

        // Insert my own
        label nTotalCells = 0;
        forAll(cWeights, celli)
        {
            allXadj[nTotalCells] = xadj[celli];
            allWeights[nTotalCells++] = cWeights[celli];
        }
        nTotalConnections = 0;
        forAll(adjncy, i)
        {
            allAdjncy[nTotalConnections++] = adjncy[i];
        }

        for (const int slave : Pstream::subProcs())
        {
            IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
            List<label> nbrAdjncy(fromSlave);
            List<label> nbrXadj(fromSlave);
            List<scalar> nbrWeights(fromSlave);

            // Append.
            forAll(nbrXadj, celli)
            {
                allXadj[nTotalCells] = nTotalConnections+nbrXadj[celli];
                allWeights[nTotalCells++] = nbrWeights[celli];
            }
            // No need to renumber xadj since already global.
            forAll(nbrAdjncy, i)
            {
                allAdjncy[nTotalConnections++] = nbrAdjncy[i];
            }
        }
        allXadj[nTotalCells] = nTotalConnections;

        labelList allDecomp;
        decomposeSerial
        (
            allAdjncy,
            allXadj,
            allWeights,
            allDecomp
        );


        // Send allFinalDecomp back
        for (const int slave : Pstream::subProcs())
        {
            OPstream toSlave(Pstream::commsTypes::scheduled, slave);
            toSlave << SubList<label>
            (
                allDecomp,
                globalCells.localSize(slave),
                globalCells.offset(slave)
            );
        }

        // Get my own part (always first)
        decomp = SubList<label>(allDecomp, globalCells.localSize());
    }
    else
    {
        // Send my part of the graph (already in global numbering)
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );
            toMaster
                << adjncy
                << SubList<label>(xadj, xadj.size()-1)
                << cWeights;
        }

        // Receive back decomposition
        IPstream fromMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );
        fromMaster >> decomp;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisLikeDecomp::metisLikeDecomp
(
    const word& derivedType,
    const dictionary& decompDict,
    const word& regionName,
    int select
)
:
    decompositionMethod(decompDict, regionName),
    coeffsDict_(findCoeffsDict(derivedType + "Coeffs", select))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::metisLikeDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
) const
{
    if (points.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Can only use this decomposition method for entire mesh" << nl
            << "and supply one coordinate (cellCentre) for every cell." << nl
            << "The number of coordinates " << points.size() << nl
            << "The number of cells in the mesh " << mesh.nCells() << nl
            << exit(FatalError);
    }

    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        identity(mesh.nCells()),
        mesh.nCells(),
        true,
        cellCells
    );

    // Decompose using default weights
    labelList decomp;
    decomposeGeneral(cellCells.m(), cellCells.offsets(), pointWeights, decomp);

    return decomp;
}


Foam::labelList Foam::metisLikeDecomp::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& agglomWeights
) const
{
    if (agglom.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Size of cell-to-coarse map " << agglom.size()
            << " differs from number of cells in mesh " << mesh.nCells()
            << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        agglom,
        agglomPoints.size(),
        true,
        cellCells
    );

    // Decompose using default weights
    labelList decomp;
    decomposeGeneral(cellCells.m(), cellCells.offsets(), agglomWeights, decomp);


    // Rework back into decomposition for original mesh
    labelList fineDistribution(agglom.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = decomp[agglom[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::metisLikeDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cellWeights
) const
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorInFunction
            << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ")." << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    CompactListList<label> cellCells(globalCellCells);

    // Decompose using default weights
    labelList decomp;
    decomposeGeneral(cellCells.m(), cellCells.offsets(), cellWeights, decomp);

    return decomp;
}

// ************************************************************************* //

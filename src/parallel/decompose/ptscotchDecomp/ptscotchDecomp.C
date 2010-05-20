/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

    From scotch forum:

    By: Francois PELLEGRINI RE: Graph mapping 'strategy' string [ reply ]
    2008-08-22 10:09 Strategy handling in Scotch is a bit tricky. In order
    not to be confused, you must have a clear view of how they are built.
    Here are some rules:

    1- Strategies are made up of "methods" which are combined by means of
    "operators".

    2- A method is of the form "m{param=value,param=value,...}", where "m"
    is a single character (this is your first error: "f" is a method name,
    not a parameter name).

    3- There exist different sort of strategies : bipartitioning strategies,
    mapping strategies, ordering strategies, which cannot be mixed. For
    instance, you cannot build a bipartitioning strategy and feed it to a
    mapping method (this is your second error).

    To use the "mapCompute" routine, you must create a mapping strategy, not
    a bipartitioning one, and so use stratGraphMap() and not
    stratGraphBipart(). Your mapping strategy should however be based on the
    "recursive bipartitioning" method ("b"). For instance, a simple (and
    hence not very efficient) mapping strategy can be :

    "b{sep=f}"

    which computes mappings with the recursive bipartitioning method "b",
    this latter using the Fiduccia-Mattheyses method "f" to compute its
    separators.

    If you want an exact partition (see your previous post), try
    "b{sep=fx}".

    However, these strategies are not the most efficient, as they do not
    make use of the multi-level framework.

    To use the multi-level framework, try for instance:

    "b{sep=m{vert=100,low=h,asc=f}x}"

    The current default mapping strategy in Scotch can be seen by using the
    "-vs" option of program gmap. It is, to date:

    b
    {
        job=t,
        map=t,
        poli=S,
        sep=
        (
            m
            {
                asc=b
                {
                    bnd=d{pass=40,dif=1,rem=1}f{move=80,pass=-1,bal=0.005},
                    org=f{move=80,pass=-1,bal=0.005},
                    width=3
                },
                low=h{pass=10}f{move=80,pass=-1,bal=0.0005},
                type=h,
                vert=80,
                rat=0.8
            }
          | m
            {
                asc=b
                {
                    bnd=d{pass=40,dif=1,rem=1}f{move=80,pass=-1,bal=0.005},
                    org=f{move=80,pass=-1,bal=0.005},
                    width=3
                },
                low=h{pass=10}f{move=80,pass=-1,bal=0.0005},
                type=h,
                vert=80,
                rat=0.8
            }
        )
    }


\*---------------------------------------------------------------------------*/

#include "ptscotchDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "OFstream.H"

extern "C"
{
#include <stdio.h>
#include "mpi.h"
#include "ptscotch.h"
}


// Hack: scotch generates floating point errors so need to switch of error
//       trapping!
#if defined(linux) || defined(linuxAMD64) || defined(linuxIA64)
#    define LINUX
#endif

#if defined(LINUX) && defined(__GNUC__)
#    define LINUX_GNUC
#endif

#ifdef LINUX_GNUC
#   ifndef __USE_GNU
#       define __USE_GNU
#   endif
#   include <fenv.h>
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ptscotchDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        ptscotchDecomp,
        dictionaryMesh
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ptscotchDecomp::check(const int retVal, const char* str)
{
    if (retVal)
    {
        FatalErrorIn("ptscotchDecomp::decompose(..)")
            << "Call to scotch routine " << str << " failed."
            << exit(FatalError);
    }
}


// Call scotch with options from dictionary.
Foam::label Foam::ptscotchDecomp::decompose
(
    List<int>& adjncy,
    List<int>& xadj,
    const scalarField& cWeights,

    List<int>& finalDecomp
)
{
//    // Dump graph
//    if (decompositionDict_.found("ptscotchCoeffs"))
//    {
//        const dictionary& scotchCoeffs =
//            decompositionDict_.subDict("ptscotchCoeffs");
//
//        if (scotchCoeffs.lookupOrDefault("writeGraph", false))
//        {
//            OFstream str(mesh_.time().path() / mesh_.name() + ".grf");
//
//            Info<< "Dumping Scotch graph file to " << str.name() << endl
//                << "Use this in combination with gpart." << endl;
//
//            label version = 0;
//            str << version << nl;
//            // Numer of vertices
//            str << xadj.size()-1 << ' ' << adjncy.size() << nl;
//            // Numbering starts from 0
//            label baseval = 0;
//            // Has weights?
//            label hasEdgeWeights = 0;
//            label hasVertexWeights = 0;
//            label numericflag = 10*hasEdgeWeights+hasVertexWeights;
//            str << baseval << ' ' << numericflag << nl;
//            for (label cellI = 0; cellI < xadj.size()-1; cellI++)
//            {
//                label start = xadj[cellI];
//                label end = xadj[cellI+1];
//                str << end-start;
//
//                for (label i = start; i < end; i++)
//                {
//                    str << ' ' << adjncy[i];
//                }
//                str << nl;
//            }
//        }
//    }

    // Strategy
    // ~~~~~~~~

    // Default.
    SCOTCH_Strat stradat;
    check(SCOTCH_stratInit(&stradat), "SCOTCH_stratInit");

    if (decompositionDict_.found("scotchCoeffs"))
    {
        const dictionary& scotchCoeffs =
            decompositionDict_.subDict("scotchCoeffs");


        string strategy;
        if (scotchCoeffs.readIfPresent("strategy", strategy))
        {
            if (debug)
            {
                Info<< "ptscotchDecomp : Using strategy " << strategy << endl;
            }
            SCOTCH_stratDgraphMap(&stradat, strategy.c_str());
            //fprintf(stdout, "S\tStrat=");
            //SCOTCH_stratSave(&stradat, stdout);
            //fprintf(stdout, "\n");
        }
    }


    // Graph
    // ~~~~~

    List<int> velotab;


    // Check for externally provided cellweights and if so initialise weights
    scalar minWeights = gMin(cWeights);
    if (cWeights.size() > 0)
    {
        if (minWeights <= 0)
        {
            WarningIn
            (
                "ptscotchDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != xadj.size()-1)
        {
            FatalErrorIn
            (
                "ptscotchDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << xadj.size()-1
                << exit(FatalError);
        }

        // Convert to integers.
        velotab.setSize(cWeights.size());
        forAll(velotab, i)
        {
            velotab[i] = int(cWeights[i]/minWeights);
        }
    }



    SCOTCH_Dgraph grafdat;
    check(SCOTCH_dgraphInit(&grafdat, MPI_COMM_WORLD), "SCOTCH_dgraphInit");
    check
    (
        SCOTCH_dgraphBuild
        (
            &grafdat,               // grafdat
            0,                      // baseval, c-style numbering
            xadj.size()-1,          // vertlocnbr, nCells
            xadj.size()-1,          // vertlocmax
            const_cast<SCOTCH_Num*>(xadj.begin()),  // vertloctab, start index per cell into
                                    // adjncy
            &xadj[1],               // vendloctab, end index  ,,

            velotab.begin(),        // veloloctab, vertex weights
            NULL,                   // vlblloctab

            adjncy.size(),          // edgelocnbr, number of arcs
            adjncy.size(),          // edgelocsiz
            adjncy.begin(),         // edgeloctab
            NULL,                   // edgegsttab
            NULL                    // edlotab, edge weights
        ),
        "SCOTCH_dgraphBuild"
    );
    check(SCOTCH_dgraphCheck(&grafdat), "SCOTCH_dgraphCheck");


    // Architecture
    // ~~~~~~~~~~~~
    // (fully connected network topology since using switch)

    SCOTCH_Arch archdat;
    check(SCOTCH_archInit(&archdat), "SCOTCH_archInit");

    List<label> processorWeights;
    if (decompositionDict_.found("scotchCoeffs"))
    {
        const dictionary& scotchCoeffs =
            decompositionDict_.subDict("scotchCoeffs");

        scotchCoeffs.readIfPresent("processorWeights", processorWeights);
    }
    if (processorWeights.size())
    {
        if (debug)
        {
            Info<< "ptscotchDecomp : Using procesor weights " << processorWeights
                << endl;
        }
        check
        (
            SCOTCH_archCmpltw(&archdat, nProcessors_, processorWeights.begin()),
            "SCOTCH_archCmpltw"
        );
    }
    else
    {
        check
        (
            SCOTCH_archCmplt(&archdat, nProcessors_),
            "SCOTCH_archCmplt"
        );
    }


    //SCOTCH_Mapping mapdat;
    //SCOTCH_dgraphMapInit(&grafdat, &mapdat, &archdat, NULL);
    //SCOTCH_dgraphMapCompute(&grafdat, &mapdat, &stradat); /*Perform mapping*/
    //SCOTCHdgraphMapExit(&grafdat, &mapdat);


    // Hack:switch off fpu error trapping
#   ifdef LINUX_GNUC
    int oldExcepts = fedisableexcept
    (
        FE_DIVBYZERO
      | FE_INVALID
      | FE_OVERFLOW
    );
#   endif

    finalDecomp.setSize(xadj.size()-1);
    finalDecomp = 0;
    check
    (
        SCOTCH_dgraphMap
        (
            &grafdat,
            &archdat,
            &stradat,           // const SCOTCH_Strat *
            finalDecomp.begin() // parttab
        ),
        "SCOTCH_graphMap"
    );

#   ifdef LINUX_GNUC
    feenableexcept(oldExcepts);
#   endif



    //finalDecomp.setSize(xadj.size()-1);
    //check
    //(
    //    SCOTCH_dgraphPart
    //    (
    //        &grafdat,
    //        nProcessors_,       // partnbr
    //        &stradat,           // const SCOTCH_Strat *
    //        finalDecomp.begin() // parttab
    //    ),
    //    "SCOTCH_graphPart"
    //);

    // Release storage for graph
    SCOTCH_dgraphExit(&grafdat);
    // Release storage for strategy
    SCOTCH_stratExit(&stradat);
    // Release storage for network topology
    SCOTCH_archExit(&archdat);

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ptscotchDecomp::ptscotchDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::ptscotchDecomp::decompose
(
    const pointField& points,
    const scalarField& pointWeights
)
{
    if (points.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "ptscotchDecomp::decompose(const pointField&, const scalarField&)"
        )
            << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << points.size() << endl
            << "The number of cells in the mesh " << mesh_.nCells()
            << exit(FatalError);
    }

//    // For running sequential ...
//    if (Pstream::nProcs() <= 1)
//    {
//        return scotchDecomp(decompositionDict_, mesh_)
//            .decompose(points, pointWeights);
//    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    // Connections
    Field<int> adjncy;
    // Offsets into adjncy
    Field<int> xadj;
    calcDistributedCSR
    (
        mesh_,
        adjncy,
        xadj
    );

    // Decompose using default weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, pointWeights, finalDecomp);

    // Copy back to labelList
    labelList decomp(finalDecomp.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


Foam::labelList Foam::ptscotchDecomp::decompose
(
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& pointWeights
)
{
    if (agglom.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "ptscotchDecomp::decompose(const labelList&, const pointField&)"
        )   << "Size of cell-to-coarse map " << agglom.size()
            << " differs from number of cells in mesh " << mesh_.nCells()
            << exit(FatalError);
    }

//    // For running sequential ...
//    if (Pstream::nProcs() <= 1)
//    {
//        return scotchDecomp(decompositionDict_, mesh_)
//            .decompose(agglom, agglomPoints, pointWeights);
//    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    List<int> adjncy;
    List<int> xadj;
    {
        // Get cellCells on coarse mesh.
        labelListList cellCells;
        calcCellCells
        (
            mesh_,
            agglom,
            agglomPoints.size(),
            cellCells
        );

        calcCSR(cellCells, adjncy, xadj);
    }

    // Decompose using weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, pointWeights, finalDecomp);

    // Rework back into decomposition for original mesh_
    labelList fineDistribution(agglom.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = finalDecomp[agglom[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::ptscotchDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
)
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorIn
        (
            "ptscotchDecomp::decompose(const pointField&, const labelListList&)"
        )   << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ")." << exit(FatalError);
    }

//    // For running sequential ...
//    if (Pstream::nProcs() <= 1)
//    {
//        return scotchDecomp(decompositionDict_, mesh_)
//            .decompose(globalCellCells, cellCentres, cWeights);
//    }


    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    List<int> adjncy;
    List<int> xadj;
    calcCSR(globalCellCells, adjncy, xadj);

    // Decompose using weights
    List<int> finalDecomp;
    decompose(adjncy, xadj, cWeights, finalDecomp);

    // Copy back to labelList
    labelList decomp(finalDecomp.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


// ************************************************************************* //

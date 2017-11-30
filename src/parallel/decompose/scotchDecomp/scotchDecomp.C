/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

    m
    {
        asc=b
        {
            width=3,
            bnd=d{pass=40, dif=1, rem=0}
            f{move=80, pass=-1, bal=0.01},
            org=f{move=80,pass=-1,bal=0.01}
        },
        low=r
        {
            job=t,
            bal=0.01,
            map=t,
            poli=S,
            sep=
            (
                m
                {
                    asc=b
                    {
                        bnd=f{move=120, pass=-1, bal=0.01, type=b},
                        org=f{move=120,pass=-1,bal=0.01,type=b},
                        width=3
                    },
                    low=h{pass=10}
                    f{move=120,pass=-1,bal=0.01,type=b},
                    vert=120,
                    rat=0.8
                }
               |m
                {
                    asc=b
                    {
                        bnd=f{move=120,pass=-1,bal=0.01,type=b},
                        org=f{move=120,pass=-1,bal=0.01,type=b},
                        width=3
                    },
                    low=h{pass=10}
                    f{move=120,pass=-1,bal=0.01,type=b},
                    vert=120,
                    rat=0.8
                }
            )
        },
        vert=10000,
        rat=0.8,
        type=0
    }


    Note: instead of gmap run gpart \<nProcs\> -vs \<grfFile\>
    where \<grfFile\> can be obtained by running with 'writeGraph=true'

\*---------------------------------------------------------------------------*/

#include "scotchDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "Time.H"
#include "OFstream.H"

extern "C"
{
#include "scotch.h"
}


// Hack: scotch generates floating point errors so need to switch of error
//       trapping!
#ifdef __GLIBC__
    #ifndef _GNU_SOURCE
        #define _GNU_SOURCE
    #endif
    #include <fenv.h>
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(scotchDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotchDecomp,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotchDecomp,
        dictionaryRegion
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::scotchDecomp::graphPath(const polyMesh& mesh)
{
    graphPath_ = mesh.time().path()/mesh.name() + ".grf";
}


void Foam::scotchDecomp::check(const int retVal, const char* str)
{
    if (retVal)
    {
        FatalErrorInFunction
            << "Call to scotch routine " << str << " failed."
            << exit(FatalError);
    }
}


Foam::label Foam::scotchDecomp::decomposeSerial
(
    const labelUList& adjncy,
    const labelUList& xadj,
    const UList<scalar>& cWeights,
    List<label>& decomp
)
{
    // Dump graph
    if (coeffsDict_.lookupOrDefault("writeGraph", false))
    {
        OFstream str(graphPath_);

        Info<< "Dumping Scotch graph file to " << str.name() << endl
            << "Use this in combination with gpart." << endl;

        const label version = 0;
        str << version << nl;
        // Numer of vertices
        str << xadj.size()-1 << ' ' << adjncy.size() << nl;

        // Numbering starts from 0
        const label baseval = 0;
        // Has weights?
        const label hasEdgeWeights = 0;
        const label hasVertexWeights = 0;
        const label numericflag = 10*hasEdgeWeights+hasVertexWeights;
        str << baseval << ' ' << numericflag << nl;

        for (label celli = 0; celli < xadj.size()-1; ++celli)
        {
            const label start = xadj[celli];
            const label end = xadj[celli+1];

            str << end-start; // size

            for (label i = start; i < end; ++i)
            {
                str << ' ' << adjncy[i];
            }
            str << nl;
        }
    }

    // Strategy
    // ~~~~~~~~

    // Default.
    SCOTCH_Strat stradat;
    check(SCOTCH_stratInit(&stradat), "SCOTCH_stratInit");

    string strategy;
    if (coeffsDict_.readIfPresent("strategy", strategy))
    {
        if (debug)
        {
            Info<< "scotchDecomp : Using strategy " << strategy << endl;
        }
        SCOTCH_stratGraphMap(&stradat, strategy.c_str());
        //fprintf(stdout, "S\tStrat=");
        //SCOTCH_stratSave(&stradat, stdout);
        //fprintf(stdout, "\n");
    }


    // Graph
    // ~~~~~

    List<label> velotab;

    // Check for externally provided cellweights and if so initialise weights
    // Note: min, not gMin since routine runs on master only.
    const scalar minWeights = min(cWeights);
    if (!cWeights.empty())
    {
        if (minWeights <= 0)
        {
            WarningInFunction
                << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != xadj.size()-1)
        {
            FatalErrorInFunction
                << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << xadj.size()-1
                << exit(FatalError);
        }

        scalar velotabSum = sum(cWeights)/minWeights;

        scalar rangeScale(1.0);

        if (velotabSum > scalar(labelMax - 1))
        {
            // 0.9 factor of safety to avoid floating point round-off in
            // rangeScale tipping the subsequent sum over the integer limit.
            rangeScale = 0.9*scalar(labelMax - 1)/velotabSum;

            WarningInFunction
                << "Sum of weights has overflowed integer: " << velotabSum
                << ", compressing weight scale by a factor of " << rangeScale
                << endl;
        }

        // Convert to integers.
        velotab.setSize(cWeights.size());

        forAll(velotab, i)
        {
            velotab[i] = int((cWeights[i]/minWeights - 1)*rangeScale) + 1;
        }
    }


    SCOTCH_Graph grafdat;
    check(SCOTCH_graphInit(&grafdat), "SCOTCH_graphInit");
    check
    (
        SCOTCH_graphBuild
        (
            &grafdat,
            0,                      // baseval, c-style numbering
            xadj.size()-1,          // vertnbr, nCells
            xadj.begin(),           // verttab, start index per cell into adjncy
            &xadj[1],               // vendtab, end index  ,,
            velotab.begin(),        // velotab, vertex weights
            nullptr,                   // vlbltab
            adjncy.size(),          // edgenbr, number of arcs
            adjncy.begin(),         // edgetab
            nullptr                    // edlotab, edge weights
        ),
        "SCOTCH_graphBuild"
    );
    check(SCOTCH_graphCheck(&grafdat), "SCOTCH_graphCheck");


    // Architecture
    // ~~~~~~~~~~~~
    // (fully connected network topology since using switch)

    SCOTCH_Arch archdat;
    check(SCOTCH_archInit(&archdat), "SCOTCH_archInit");

    List<label> processorWeights;
    if
    (
        coeffsDict_.readIfPresent("processorWeights", processorWeights)
     && processorWeights.size()
    )
    {
        if (debug)
        {
            Info<< "scotchDecomp : Using procesor weights " << processorWeights
                << endl;
        }
        check
        (
            SCOTCH_archCmpltw
            (
                &archdat, nDomains_, processorWeights.begin()
            ),
            "SCOTCH_archCmpltw"
        );
    }
    else
    {
        check
        (
            SCOTCH_archCmplt(&archdat, nDomains_),
            "SCOTCH_archCmplt"
        );


        //- Hack to test clustering. Note that decomp is non-compact
        //  numbers!
        //
        ////- Set up variable sizes architecture
        //check
        //(
        //    SCOTCH_archVcmplt(&archdat),
        //    "SCOTCH_archVcmplt"
        //);
        //
        ////- Stategy flags: go for quality or load balance (or leave default)
        //SCOTCH_Num straval = 0;
        ////straval |= SCOTCH_STRATQUALITY;
        ////straval |= SCOTCH_STRATQUALITY;
        //
        ////- Number of cells per agglomeration
        ////SCOTCH_Num agglomSize = SCOTCH_archSize(&archdat);
        //SCOTCH_Num agglomSize = 3;
        //
        ////- Build strategy for agglomeration
        //check
        //(
        //    SCOTCH_stratGraphClusterBuild
        //    (
        //        &stradat,   // strategy to build
        //        straval,    // strategy flags
        //        agglomSize, // cells per cluster
        //        1.0,        // weight?
        //        0.01        // max load imbalance
        //    ),
        //    "SCOTCH_stratGraphClusterBuild"
        //);
    }


    //SCOTCH_Mapping mapdat;
    //SCOTCH_graphMapInit(&grafdat, &mapdat, &archdat, nullptr);
    //SCOTCH_graphMapCompute(&grafdat, &mapdat, &stradat); /* Perform mapping */
    //SCOTCH_graphMapExit(&grafdat, &mapdat);


    // Hack:switch off fpu error trapping
    #ifdef FE_NOMASK_ENV
    int oldExcepts = fedisableexcept
    (
        FE_DIVBYZERO
      | FE_INVALID
      | FE_OVERFLOW
    );
    #endif

    decomp.setSize(xadj.size()-1);
    decomp = 0;
    check
    (
        SCOTCH_graphMap
        (
            &grafdat,
            &archdat,
            &stradat,       // const SCOTCH_Strat *
            decomp.begin()  // parttab
        ),
        "SCOTCH_graphMap"
    );

    #ifdef FE_NOMASK_ENV
    feenableexcept(oldExcepts);
    #endif

    //decomp.setSize(xadj.size()-1);
    //check
    //(
    //    SCOTCH_graphPart
    //    (
    //        &grafdat,
    //        nDomains_,      // partnbr
    //        &stradat,       // const SCOTCH_Strat *
    //        decomp.begin()  // parttab
    //    ),
    //    "SCOTCH_graphPart"
    //);

    // Release storage for graph
    SCOTCH_graphExit(&grafdat);
    // Release storage for strategy
    SCOTCH_stratExit(&stradat);
    // Release storage for network topology
    SCOTCH_archExit(&archdat);

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scotchDecomp::scotchDecomp(const dictionary& decompDict)
:
    metisLikeDecomp(typeName, decompDict, selectionType::NULL_DICT)
{}


Foam::scotchDecomp::scotchDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    metisLikeDecomp(typeName, decompDict, regionName, selectionType::NULL_DICT)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::scotchDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{
    // Where to write graph
    graphPath(mesh);

    return metisLikeDecomp::decompose
    (
        mesh,
        points,
        pointWeights
    );
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& pointWeights
)
{
    // Where to write graph
    graphPath(mesh);

    return metisLikeDecomp::decompose
    (
        mesh,
        agglom,
        agglomPoints,
        pointWeights
    );
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
)
{
    // Where to write graph
    graphPath_ = "scotch.grf";

    return metisLikeDecomp::decompose
    (
        globalCellCells,
        cellCentres,
        cWeights
    );
}


// ************************************************************************* //

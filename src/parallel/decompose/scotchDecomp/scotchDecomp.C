/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "scotchDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "Time.H"
#include "OFstream.H"

// Probably not needed, but in case we pickup a ptscotch.h ...
#define MPICH_SKIP_MPICXX
#define OMPI_SKIP_MPICXX

#include "scotch.h"

// Hack: scotch generates floating point errors so need to switch off error
//       trapping!
#ifdef __GLIBC__
    #ifndef _GNU_SOURCE
        #define _GNU_SOURCE
    #endif
    #include <fenv.h>
#endif

// Provide a clear error message if we have a size mismatch
static_assert
(
    sizeof(Foam::label) == sizeof(SCOTCH_Num),
    "sizeof(Foam::label) == sizeof(SCOTCH_Num), check your scotch headers"
);


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(scotchDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotchDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::scotchDecomp::graphPath(const polyMesh& mesh) const
{
    graphPath_ = mesh.time().path()/mesh.name() + ".grf";
}


void Foam::scotchDecomp::check(const int retVal, const char* str)
{
    if (retVal)
    {
        FatalErrorInFunction
            << "Call to scotch routine " << str << " failed.\n"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::scotchDecomp::decomposeSerial
(
    const labelList& adjncy,
    const labelList& xadj,
    const List<scalar>& cWeights,
    labelList& decomp
) const
{
    // Dump graph
    if (coeffsDict_.getOrDefault("writeGraph", false))
    {
        OFstream str(graphPath_);

        Info<< "Dumping Scotch graph file to " << str.name() << nl
            << "Use this in combination with gpart." << endl;

        const label version = 0;
        str << version << nl;
        // Number of vertices
        str << xadj.size()-1 << ' ' << adjncy.size() << nl;

        // Numbering starts from 0
        const label baseval = 0;

        // Has weights?
        const label hasEdgeWeights = 0;
        const label hasVertexWeights = 0;
        const label numericflag = 10*hasEdgeWeights+hasVertexWeights;
        str << baseval << ' ' << numericflag << nl;

        for (label celli = 1; celli < xadj.size(); ++celli)
        {
            const label start = xadj[celli-1];
            const label end = xadj[celli];

            str << end-start; // size

            for (label i = start; i < end; ++i)
            {
                str << ' ' << adjncy[i];
            }
            str << nl;
        }
    }

    // Make repeatable
    SCOTCH_randomReset();

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

    labelList velotab;

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

    labelList processorWeights;
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
        if (processorWeights.size() != nDomains_)
        {
            FatalIOErrorInFunction(coeffsDict_)
                << "processorWeights not the same size"
                << " as the wanted number of domains " << nDomains_
                << exit(FatalIOError);
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
) const
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
) const
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
) const
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

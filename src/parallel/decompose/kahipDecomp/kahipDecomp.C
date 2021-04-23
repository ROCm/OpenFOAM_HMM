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

#include "kahipDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "PrecisionAdaptor.H"

#include "kaHIP_interface.h"

#include <string>
#include <map>
#include <vector>

// Provide a clear error message if we have a severe size mismatch
// Allow widening, but not narrowing

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kahipDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        kahipDecomp,
        dictionary
    );
}


const Foam::Enum
<
    Foam::kahipDecomp::configs
>
Foam::kahipDecomp::configNames
({
    { kahipDecomp::configs::FAST, "fast" },
    { kahipDecomp::configs::ECO, "eco" },
    { kahipDecomp::configs::STRONG, "strong" },
    { kahipDecomp::configs::FASTSOCIAL, "fast-social" },
    { kahipDecomp::configs::ECOSOCIAL, "eco-social" },
    { kahipDecomp::configs::STRONGSOCIAL, "strong-social" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::kahipDecomp::decomposeSerial
(
    const labelList& adjncy,
    const labelList& xadj,
    const List<scalar>& cWeights,
    labelList& decomp
) const
{
    // Default setup
    enum configs kahipConfig = configs::FAST;
    double imbalance = 0.01;
    int seed = 0;
    bool verbose = false;

    #if WM_LABEL_SIZE == 64
    if (xadj.size()-1 > INT_MAX)
    {
        FatalErrorInFunction
            << "Cannot decompose " << (xadj.size()-1) << " cells," << nl
            << "Exceeded integer limit of " << INT_MAX << nl
            << exit(FatalError);
    }
    #endif

    int numCells = xadj.size()-1;

    // Cell weights (so on the vertices of the dual)
    List<int> cellWeights;

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

        if (cWeights.size() != numCells)
        {
            FatalErrorInFunction
                << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << numCells
                << exit(FatalError);
        }

        // Convert to integers.
        cellWeights.setSize(cWeights.size());
        forAll(cellWeights, i)
        {
            cellWeights[i] = int(cWeights[i]/minWeights);
        }
    }

    configNames.readIfPresent("config", coeffsDict_, kahipConfig);
    coeffsDict_.readIfPresent("imbalance", imbalance);
    coeffsDict_.readIfPresent("verbose", verbose);

    Info<< "kahipDecomp :"
        << " config=" << configNames[kahipConfig]
        << " imbalance=" << imbalance;

    if (coeffsDict_.readIfPresent("seed", seed))
    {
        Info<< " seed=" << seed;
    }

    // Additional sizing parameters (testing only)
    std::map<std::string, std::vector<int>> sizingParams;

    List<int> labels;
    if
    (
        coeffsDict_.readIfPresent("hierarchy", labels)
     && !labels.empty()
    )
    {
        std::vector<int> vec;
        vec.reserve(labels.size()+1);

        // Verify sizing

        int n = 1;
        for (const auto val : labels)
        {
            n *= val;
            vec.push_back(val);
        }

        if (n != nDomains_)
        {
            // Size mismatch. Try to correct.

            if (nDomains_ % n)
            {
                WarningInFunction
                    << "Mismatch in number of processors and "
                    << "hierarchy specified" << flatOutput(labels) << endl;

                vec.clear();
            }
            else
            {
                // Evenly divisible, add extra hierarchy level
                vec.push_back(nDomains_ / n);
            }
        }

        if (!vec.empty())
        {
            sizingParams["hierarchy"] = std::move(vec);
            Info<< " hierarchy=" << flatOutput(labels);
        }
    }

    if
    (
        coeffsDict_.readIfPresent("distance", labels)
     && !labels.empty()
    )
    {
        std::vector<int> vec(labels.size());

        forAll(labels, i)
        {
            vec[i] = labels[i];
        }

        sizingParams["distance"] = std::move(vec);
        Info<< " distance=" << flatOutput(labels);
    }

    Info<< endl;


    // Number of partitions
    int nParts = nDomains_;

    // Output: number of cut edges
    int edgeCut = 0;

    // Addressing
    ConstPrecisionAdaptor<int, label, List> xadj_param(xadj);
    ConstPrecisionAdaptor<int, label, List> adjncy_param(adjncy);

    // Output: cell -> processor addressing
    decomp.resize(numCells);
    PrecisionAdaptor<int, label, List> decomp_param(decomp, false);


#if 0 // WIP: #ifdef KAFFPA_CPP_INTERFACE
    kaffpa_cpp
    (
        &numCells,          // num vertices in graph
        (cellWeights.size() ? cellWeights.begin() : nullptr), // vertex wts
        xadj_param.constCast().data(),          // indexing into adjncy
        nullptr,            // edge wts
        adjncy_param.constCast().data(),        // neighbour info
        &nParts,            // nparts
        &imbalance,         // amount of imbalance allowed
        !verbose,           // suppress output
        seed,               // for random
        int(kahipConfig),
        &edgeCut,                   // [output]
        decomp_param.ref().data(),  // [output]
        sizingParams
    );
#else
    kaffpa
    (
        &numCells,          // num vertices in graph
        (cellWeights.size() ? cellWeights.begin() : nullptr), // vertex wts
        xadj_param.constCast().data(),          // indexing into adjncy
        nullptr,            // edge wts
        adjncy_param.constCast().data(),        // neighbour info
        &nParts,            // nparts
        &imbalance,         // amount of imbalance allowed
        !verbose,           // suppress output
        seed,               // for random
        int(kahipConfig),
        &edgeCut,                   // [output]
        decomp_param.ref().data()   // [output]
    );
#endif

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kahipDecomp::kahipDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    metisLikeDecomp(typeName, decompDict, regionName, selectionType::NULL_DICT)
{}


// ************************************************************************* //

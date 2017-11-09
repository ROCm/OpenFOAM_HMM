/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "metisDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

extern "C"
{
    #define OMPI_SKIP_MPICXX
    #include "metis.h"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(metisDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        metisDecomp,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        metisDecomp,
        dictionaryRegion
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::metisDecomp::decomposeSerial
(
    const labelUList& adjncy,
    const labelUList& xadj,
    const UList<scalar>& cWeights,
    List<label>& decomp
)
{
    // Method of decomposition
    // recursive: multi-level recursive bisection (default)
    // k-way: multi-level k-way
    word method("recursive");

    const dictionary* coeffsDictPtr =
        decompositionDict_.subDictPtr("metisCoeffs");

    label numCells = xadj.size()-1;

    // Decomposition options
    List<label> options(METIS_NOPTIONS);
    METIS_SetDefaultOptions(options.begin());

    // Processor weights initialised with no size, only used if specified in
    // a file
    Field<real_t> processorWeights;

    // Cell weights (so on the vertices of the dual)
    List<label> cellWeights;

    // Face weights (so on the edges of the dual)
    List<label> faceWeights;

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


    // Check for user supplied weights and decomp options
    if (coeffsDictPtr)
    {
        const dictionary& coeffDict = *coeffsDictPtr;

        word weightsFile;

        if (coeffDict.readIfPresent("method", method))
        {
            if (method != "recursive" && method != "k-way")
            {
                FatalErrorInFunction
                    << "Method " << method << " in metisCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be 'recursive' or 'k-way'"
                    << exit(FatalError);
            }

            Info<< "metisDecomp : Using Metis method     " << method
                << nl << endl;
        }

        if (coeffDict.readIfPresent("options", options))
        {
            if (options.size() != METIS_NOPTIONS)
            {
                FatalErrorInFunction
                    << "Number of options in metisCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be " << METIS_NOPTIONS
                    << exit(FatalError);
            }

            Info<< "metisDecomp : Using Metis options     " << options
                << nl << endl;
        }

        if (coeffDict.readIfPresent("processorWeights", processorWeights))
        {
            processorWeights /= sum(processorWeights);

            if (processorWeights.size() != nDomains_)
            {
                FatalErrorInFunction
                    << "Number of processor weights "
                    << processorWeights.size()
                    << " does not equal number of domains " << nDomains_
                    << exit(FatalError);
            }
        }
    }

    label ncon = 1;
    label nProcs = nDomains_;

    // Output: cell -> processor addressing
    decomp.setSize(numCells);

    // Output: number of cut edges
    label edgeCut = 0;

    if (method == "recursive")
    {
        METIS_PartGraphRecursive
        (
            &numCells,          // num vertices in graph
            &ncon,              // num balancing constraints
            const_cast<labelUList&>(xadj).begin(),   // indexing into adjncy
            const_cast<labelUList&>(adjncy).begin(), // neighbour info
            cellWeights.begin(),// vertex wts
            nullptr,               // vsize: total communication vol
            faceWeights.begin(),// edge wts
            &nProcs,            // nParts
            processorWeights.begin(),   // tpwgts
            nullptr,               // ubvec: processor imbalance (default)
            options.begin(),
            &edgeCut,
            decomp.begin()
        );
    }
    else
    {
        METIS_PartGraphKway
        (
            &numCells,          // num vertices in graph
            &ncon,              // num balancing constraints
            const_cast<labelUList&>(xadj).begin(),   // indexing into adjncy
            const_cast<labelUList&>(adjncy).begin(), // neighbour info
            cellWeights.begin(),// vertex wts
            nullptr,               // vsize: total communication vol
            faceWeights.begin(),// edge wts
            &nProcs,            // nParts
            processorWeights.begin(),   // tpwgts
            nullptr,               // ubvec: processor imbalance (default)
            options.begin(),
            &edgeCut,
            decomp.begin()
        );
    }

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisDecomp::metisDecomp(const dictionary& decompDict)
:
    metisLikeDecomp(typeName, decompDict, selectionType::NULL_DICT)
{}


Foam::metisDecomp::metisDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    metisLikeDecomp(typeName, decompDict, regionName, selectionType::NULL_DICT)
{}


// ************************************************************************* //

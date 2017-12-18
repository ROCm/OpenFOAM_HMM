/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "decompositionMethod.H"
#include "globalIndex.H"
#include "syncTools.H"
#include "Tuple2.H"
#include "faceSet.H"
#include "regionSplit.H"
#include "localPointRegion.H"
#include "minData.H"
#include "FaceCellWave.H"

#include "preserveBafflesConstraint.H"
#include "preservePatchesConstraint.H"
#include "preserveFaceZonesConstraint.H"
#include "singleProcessorFaceSetsConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decompositionMethod, 0);
    defineRunTimeSelectionTable(decompositionMethod, dictionary);
    defineRunTimeSelectionTable(decompositionMethod, dictionaryRegion);

    // Fallback name when searching for optional coefficients directories
    static const word defaultName("coeffs");
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::decompositionMethod::nDomains(const dictionary& decompDict)
{
    return readLabel(decompDict.lookup("numberOfSubdomains"));
}


Foam::label Foam::decompositionMethod::nDomains
(
    const dictionary& decompDict,
    const word& regionName
)
{
    const label nDomainsGlobal = nDomains(decompDict);

    const dictionary& regionDict(optionalRegionDict(decompDict, regionName));

    label nDomainsRegion;
    if (regionDict.readIfPresent("numberOfSubdomains", nDomainsRegion))
    {
        if (nDomainsRegion >= 1 && nDomainsRegion <= nDomainsGlobal)
        {
            return nDomainsRegion;
        }

        WarningInFunction
            << "ignoring out of range numberOfSubdomains "
            << nDomainsRegion << " for region " << regionName
            << nl << nl
            << endl;
    }

    return nDomainsGlobal;
}


const Foam::dictionary& Foam::decompositionMethod::optionalRegionDict
(
    const dictionary& decompDict,
    const word& regionName
)
{
    auto finder = decompDict.csearch("regions");

    if (!regionName.empty() && finder.isDict())
    {
        finder = finder.dict().csearch(regionName);

        if (finder.isDict())
        {
            return finder.dict();
        }
    }

    return dictionary::null;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::decompositionMethod::readConstraints()
{
    constraints_.clear();

    // Read any constraints
    wordList constraintTypes;

    const dictionary* dictptr = decompositionDict_.subDictPtr("constraints");

    if (dictptr)
    {
        forAllConstIters(*dictptr, iter)
        {
            const dictionary& dict = iter().dict();

            constraintTypes.append(dict.lookup("type"));

            constraints_.append
            (
                decompositionConstraint::New
                (
                    dict,
                    constraintTypes.last()
                )
            );
        }
    }

    // Backwards compatibility
    if
    (
        decompositionDict_.found("preserveBaffles")
     && !constraintTypes.found
        (
            decompositionConstraints::preserveBafflesConstraint::typeName
        )
    )
    {
        constraints_.append
        (
            new decompositionConstraints::preserveBafflesConstraint()
        );
    }

    if
    (
        decompositionDict_.found("preservePatches")
     && !constraintTypes.found
        (
            decompositionConstraints::preservePatchesConstraint::typeName
        )
    )
    {
        const wordReList pNames(decompositionDict_.lookup("preservePatches"));

        constraints_.append
        (
            new decompositionConstraints::preservePatchesConstraint(pNames)
        );
    }

    if
    (
        decompositionDict_.found("preserveFaceZones")
     && !constraintTypes.found
        (
            decompositionConstraints::preserveFaceZonesConstraint::typeName
        )
    )
    {
        const wordReList zNames(decompositionDict_.lookup("preserveFaceZones"));

        constraints_.append
        (
            new decompositionConstraints::preserveFaceZonesConstraint(zNames)
        );
    }

    if
    (
        decompositionDict_.found("singleProcessorFaceSets")
     && !constraintTypes.found
        (
            decompositionConstraints::preserveFaceZonesConstraint::typeName
        )
    )
    {
        const List<Tuple2<word, label>> zNameAndProcs
        (
            decompositionDict_.lookup("singleProcessorFaceSets")
        );

        constraints_.append
        (
            new decompositionConstraints::singleProcessorFaceSetsConstraint
            (
                zNameAndProcs
            )
        );
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::dictionary& Foam::decompositionMethod::findCoeffsDict
(
    const dictionary& dict,
    const word& coeffsName,
    int select
)
{
    dictionary::const_searcher fnd;

    if
    (
        (fnd = dict.csearch(coeffsName)).isDict()
     ||
        (
            !(select & selectionType::EXACT)
         && (fnd = dict.csearch(defaultName)).isDict()
        )
    )
    {
        return fnd.dict();
    }

    // Not found
    if (select & selectionType::MANDATORY)
    {
        FatalIOError
            << "'" << coeffsName << "' dictionary not found in dictionary "
            << dict.name() << endl
            << abort(FatalIOError);
    }

    if (select & selectionType::NULL_DICT)
    {
        return dictionary::null;
    }

    return dict;
}


const Foam::dictionary& Foam::decompositionMethod::findCoeffsDict
(
    const word& coeffsName,
    int select
) const
{
    dictionary::const_searcher fnd;

    if
    (
        !decompositionRegionDict_.empty()
    &&
        (
            (fnd = decompositionRegionDict_.csearch(coeffsName)).isDict()
         ||
            (
                !(select & selectionType::EXACT)
             && (fnd = decompositionRegionDict_.csearch(defaultName)).isDict()
            )
        )
    )
    {
        return fnd.dict();
    }

    if
    (
        (fnd = decompositionDict_.csearch(coeffsName)).isDict()
     ||
        (
            !(select & selectionType::EXACT)
         && (fnd = decompositionDict_.csearch(defaultName)).isDict()
        )
    )
    {
        return fnd.dict();
    }

    // Not found
    if (select & selectionType::MANDATORY)
    {
        FatalIOError
            << "'" << coeffsName << "' dictionary not found in dictionary "
            << decompositionDict_.name() << endl
            << abort(FatalIOError);
    }

    if (select & selectionType::NULL_DICT)
    {
        return dictionary::null;
    }

    return decompositionDict_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethod::decompositionMethod
(
    const dictionary& decompDict
)
:
    decompositionDict_(decompDict),
    decompositionRegionDict_(dictionary::null),
    nDomains_(nDomains(decompDict))
{
    readConstraints();
}


Foam::decompositionMethod::decompositionMethod
(
    const dictionary& decompDict,
    const word& regionName
)
:
    decompositionDict_(decompDict),
    decompositionRegionDict_
    (
        optionalRegionDict(decompositionDict_, regionName)
    ),
    nDomains_(nDomains(decompDict, regionName))
{
    readConstraints();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::decompositionMethod> Foam::decompositionMethod::New
(
    const dictionary& decompDict
)
{
    const word methodType(decompDict.lookup("method"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(methodType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown decompositionMethod "
            << methodType << nl << nl
            << "Valid decompositionMethods : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // verbose
    {
        Info<< "Selecting decompositionMethod " << methodType
            << " [" << (nDomains(decompDict)) << "]" << endl;
    }

    return autoPtr<decompositionMethod>(cstrIter()(decompDict));
}


Foam::autoPtr<Foam::decompositionMethod> Foam::decompositionMethod::New
(
    const dictionary& decompDict,
    const word& regionName

)
{
    const dictionary& regionDict(optionalRegionDict(decompDict, regionName));

    if (regionDict.empty())
    {
        // No region-specific information - just forward to normal routine
        return decompositionMethod::New(decompDict);
    }

    word methodType(decompDict.lookup("method"));
    regionDict.readIfPresent("method", methodType);

    auto cstrIter = dictionaryRegionConstructorTablePtr_->cfind(methodType);

    if (!cstrIter.found())
    {
        WarningInFunction
            << nl
            << "Unknown region decompositionMethod "
            << methodType << nl << nl
            << "Valid decompositionMethods : " << endl
            << dictionaryRegionConstructorTablePtr_->sortedToc() << nl
            << "Reverting to non-region version" << nl
            << endl;

        return decompositionMethod::New(decompDict);
    }

    // verbose
    {
        Info<< "Selecting decompositionMethod " << methodType
            << " [" << (nDomains(decompDict, regionName)) << "] (region "
            << regionName << ")" << endl;
    }

    return autoPtr<decompositionMethod>(cstrIter()(decompDict, regionName));
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
        true,                       // use global cell labels
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


void Foam::decompositionMethod::calcCellCells
(
    const polyMesh& mesh,
    const labelList& agglom,
    const label nLocalCoarse,
    const bool parallel,
    CompactListList<label>& cellCells
)
{
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Create global cell numbers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    globalIndex globalAgglom
    (
        nLocalCoarse,
        Pstream::msgType(),
        Pstream::worldComm,
        parallel
    );


    // Get agglomerate owner on other side of coupled faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList globalNeighbour(mesh.nFaces()-mesh.nInternalFaces());

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
        {
            label facei = pp.start();
            label bFacei = pp.start() - mesh.nInternalFaces();

            forAll(pp, i)
            {
                globalNeighbour[bFacei] = globalAgglom.toGlobal
                (
                    agglom[faceOwner[facei]]
                );

                bFacei++;
                facei++;
            }
        }
    }

    // Get the cell on the other side of coupled patches
    syncTools::swapBoundaryFaceList(mesh, globalNeighbour);


    // Count number of faces (internal + coupled)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Number of faces per coarse cell
    labelList nFacesPerCell(nLocalCoarse, 0);

    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label own = agglom[faceOwner[facei]];
        const label nei = agglom[faceNeighbour[facei]];

        nFacesPerCell[own]++;
        nFacesPerCell[nei]++;
    }

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
        {
            label facei = pp.start();
            label bFacei = pp.start()-mesh.nInternalFaces();

            forAll(pp, i)
            {
                const label own = agglom[faceOwner[facei]];

                label globalNei = globalNeighbour[bFacei];
                if
                (
                   !globalAgglom.isLocal(globalNei)
                 || globalAgglom.toLocal(globalNei) != own
                )
                {
                    nFacesPerCell[own]++;
                }

                facei++;
                bFacei++;
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
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label own = agglom[faceOwner[facei]];
        const label nei = agglom[faceNeighbour[facei]];

        m[offsets[own] + nFacesPerCell[own]++] = globalAgglom.toGlobal(nei);
        m[offsets[nei] + nFacesPerCell[nei]++] = globalAgglom.toGlobal(own);
    }

    // For boundary faces is offsetted coupled neighbour
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
        {
            label facei = pp.start();
            label bFacei = pp.start()-mesh.nInternalFaces();

            forAll(pp, i)
            {
                const label own = agglom[faceOwner[facei]];

                const label globalNei = globalNeighbour[bFacei];

                if
                (
                   !globalAgglom.isLocal(globalNei)
                 || globalAgglom.toLocal(globalNei) != own
                )
                {
                    m[offsets[own] + nFacesPerCell[own]++] = globalNei;
                }

                facei++;
                bFacei++;
            }
        }
    }


    // Check for duplicates connections between cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Done as postprocessing step since we now have cellCells.
    label newIndex = 0;
    labelHashSet nbrCells;


    if (cellCells.size() == 0)
    {
        return;
    }

    label startIndex = cellCells.offsets()[0];

    forAll(cellCells, celli)
    {
        nbrCells.clear();
        nbrCells.insert(globalAgglom.toGlobal(celli));

        const label endIndex = cellCells.offsets()[celli+1];

        for (label i = startIndex; i < endIndex; ++i)
        {
            if (nbrCells.insert(cellCells.m()[i]))
            {
                cellCells.m()[newIndex++] = cellCells.m()[i];
            }
        }
        startIndex = endIndex;
        cellCells.offsets()[celli+1] = newIndex;
    }

    cellCells.m().setSize(newIndex);

    //forAll(cellCells, celli)
    //{
    //    Pout<< "Original: Coarse cell " << celli << endl;
    //    forAll(mesh.cellCells()[celli], i)
    //    {
    //        Pout<< "    nbr:" << mesh.cellCells()[celli][i] << endl;
    //    }
    //    Pout<< "Compacted: Coarse cell " << celli << endl;
    //    const labelUList cCells = cellCells[celli];
    //    forAll(cCells, i)
    //    {
    //        Pout<< "    nbr:" << cCells[i] << endl;
    //    }
    //}
}


void Foam::decompositionMethod::calcCellCells
(
    const polyMesh& mesh,
    const labelList& agglom,
    const label nLocalCoarse,
    const bool parallel,
    CompactListList<label>& cellCells,
    CompactListList<scalar>& cellCellWeights
)
{
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Create global cell numbers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    globalIndex globalAgglom
    (
        nLocalCoarse,
        Pstream::msgType(),
        Pstream::worldComm,
        parallel
    );


    // Get agglomerate owner on other side of coupled faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList globalNeighbour(mesh.nFaces()-mesh.nInternalFaces());

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
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
    labelList nFacesPerCell(nLocalCoarse, 0);

    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label own = agglom[faceOwner[facei]];
        const label nei = agglom[faceNeighbour[facei]];

        nFacesPerCell[own]++;
        nFacesPerCell[nei]++;
    }

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
        {
            label faceI = pp.start();
            label bFaceI = pp.start()-mesh.nInternalFaces();

            forAll(pp, i)
            {
                const label own = agglom[faceOwner[faceI]];

                const label globalNei = globalNeighbour[bFaceI];
                if
                (
                   !globalAgglom.isLocal(globalNei)
                 || globalAgglom.toLocal(globalNei) != own
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
    cellCellWeights.setSize(nFacesPerCell);

    nFacesPerCell = 0;

    labelList& m = cellCells.m();
    scalarList& w = cellCellWeights.m();
    const labelList& offsets = cellCells.offsets();

    // For internal faces is just offsetted owner and neighbour
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label own = agglom[faceOwner[facei]];
        const label nei = agglom[faceNeighbour[facei]];

        const label ownIndex = offsets[own] + nFacesPerCell[own]++;
        const label neiIndex = offsets[nei] + nFacesPerCell[nei]++;

        m[ownIndex] = globalAgglom.toGlobal(nei);
        w[ownIndex] = mag(mesh.faceAreas()[facei]);
        m[neiIndex] = globalAgglom.toGlobal(own);
        w[ownIndex] = mag(mesh.faceAreas()[facei]);
    }

    // For boundary faces is offsetted coupled neighbour
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
        {
            label faceI = pp.start();
            label bFaceI = pp.start()-mesh.nInternalFaces();

            forAll(pp, i)
            {
                const label own = agglom[faceOwner[faceI]];

                const label globalNei = globalNeighbour[bFaceI];

                if
                (
                   !globalAgglom.isLocal(globalNei)
                 || globalAgglom.toLocal(globalNei) != own
                )
                {
                    label ownIndex = offsets[own] + nFacesPerCell[own]++;
                    m[ownIndex] = globalNei;
                    w[ownIndex] = mag(mesh.faceAreas()[faceI]);
                }

                faceI++;
                bFaceI++;
            }
        }
    }


    // Check for duplicates connections between cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Done as postprocessing step since we now have cellCells.
    label newIndex = 0;
    labelHashSet nbrCells;


    if (cellCells.size() == 0)
    {
        return;
    }

    label startIndex = cellCells.offsets()[0];

    forAll(cellCells, celli)
    {
        nbrCells.clear();
        nbrCells.insert(globalAgglom.toGlobal(celli));

        const label endIndex = cellCells.offsets()[celli+1];

        for (label i = startIndex; i < endIndex; ++i)
        {
            if (nbrCells.insert(cellCells.m()[i]))
            {
                cellCells.m()[newIndex] = cellCells.m()[i];
                cellCellWeights.m()[newIndex] = cellCellWeights.m()[i];
                newIndex++;
            }
        }
        startIndex = endIndex;
        cellCells.offsets()[celli+1] = newIndex;
        cellCellWeights.offsets()[celli+1] = newIndex;
    }

    cellCells.m().setSize(newIndex);
    cellCellWeights.m().setSize(newIndex);
}


//void Foam::decompositionMethod::calcCellCells
//(
//    const polyMesh& mesh,
//    const boolList& blockedFace,
//    const List<labelPair>& explicitConnections,
//    const labelList& agglom,
//    const label nLocalCoarse,
//    const bool parallel,
//    CompactListList<label>& cellCells
//)
//{
//    const labelList& faceOwner = mesh.faceOwner();
//    const labelList& faceNeighbour = mesh.faceNeighbour();
//    const polyBoundaryMesh& patches = mesh.boundaryMesh();
//
//
//    // Create global cell numbers
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//    globalIndex globalAgglom
//    (
//        nLocalCoarse,
//        Pstream::msgType(),
//        Pstream::worldComm,
//        parallel
//    );
//
//
//    // Get agglomerate owner on other side of coupled faces
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//    labelList globalNeighbour(mesh.nFaces()-mesh.nInternalFaces());
//
//    forAll(patches, patchi)
//    {
//        const polyPatch& pp = patches[patchi];
//
//        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
//        {
//            label facei = pp.start();
//            label bFacei = pp.start() - mesh.nInternalFaces();
//
//            forAll(pp, i)
//            {
//                globalNeighbour[bFacei] = globalAgglom.toGlobal
//                (
//                    agglom[faceOwner[facei]]
//                );
//
//                bFacei++;
//                facei++;
//            }
//        }
//    }
//
//    // Get the cell on the other side of coupled patches
//    syncTools::swapBoundaryFaceList(mesh, globalNeighbour);
//
//
//    // Count number of faces (internal + coupled)
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//    // Number of faces per coarse cell
//    labelList nFacesPerCell(nLocalCoarse, 0);
//
//    // 1. Internal faces
//    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
//    {
//        if (!blockedFace[facei])
//        {
//            label own = agglom[faceOwner[facei]];
//            label nei = agglom[faceNeighbour[facei]];
//
//            nFacesPerCell[own]++;
//            nFacesPerCell[nei]++;
//        }
//    }
//
//    // 2. Coupled faces
//    forAll(patches, patchi)
//    {
//        const polyPatch& pp = patches[patchi];
//
//        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
//        {
//            label facei = pp.start();
//            label bFacei = pp.start()-mesh.nInternalFaces();
//
//            forAll(pp, i)
//            {
//                if (!blockedFace[facei])
//                {
//                    label own = agglom[faceOwner[facei]];
//
//                    label globalNei = globalNeighbour[bFacei];
//                    if
//                    (
//                       !globalAgglom.isLocal(globalNei)
//                     || globalAgglom.toLocal(globalNei) != own
//                    )
//                    {
//                        nFacesPerCell[own]++;
//                    }
//
//                    facei++;
//                    bFacei++;
//                }
//            }
//        }
//    }
//
//    // 3. Explicit connections between non-coupled boundary faces
//    forAll(explicitConnections, i)
//    {
//        const labelPair& baffle = explicitConnections[i];
//        label f0 = baffle.first();
//        label f1 = baffle.second();
//
//        if (!blockedFace[f0] && blockedFace[f1])
//        {
//            label f0Own = agglom[faceOwner[f0]];
//            label f1Own = agglom[faceOwner[f1]];
//
//            // Always count the connection between the two owner sides
//            if (f0Own != f1Own)
//            {
//                nFacesPerCell[f0Own]++;
//                nFacesPerCell[f1Own]++;
//            }
//
//            // Add any neighbour side connections
//            if (mesh.isInternalFace(f0))
//            {
//                label f0Nei = agglom[faceNeighbour[f0]];
//
//                if (mesh.isInternalFace(f1))
//                {
//                    // Internal faces
//                    label f1Nei = agglom[faceNeighbour[f1]];
//
//                    if (f0Own != f1Nei)
//                    {
//                        nFacesPerCell[f0Own]++;
//                        nFacesPerCell[f1Nei]++;
//                    }
//                    if (f0Nei != f1Own)
//                    {
//                        nFacesPerCell[f0Nei]++;
//                        nFacesPerCell[f1Own]++;
//                    }
//                    if (f0Nei != f1Nei)
//                    {
//                        nFacesPerCell[f0Nei]++;
//                        nFacesPerCell[f1Nei]++;
//                    }
//                }
//                else
//                {
//                    // f1 boundary face
//                    if (f0Nei != f1Own)
//                    {
//                        nFacesPerCell[f0Nei]++;
//                        nFacesPerCell[f1Own]++;
//                    }
//                }
//            }
//            else
//            {
//                if (mesh.isInternalFace(f1))
//                {
//                    label f1Nei = agglom[faceNeighbour[f1]];
//                    if (f0Own != f1Nei)
//                    {
//                        nFacesPerCell[f0Own]++;
//                        nFacesPerCell[f1Nei]++;
//                    }
//                }
//            }
//        }
//    }
//
//
//    // Fill in offset and data
//    // ~~~~~~~~~~~~~~~~~~~~~~~
//
//    cellCells.setSize(nFacesPerCell);
//
//    nFacesPerCell = 0;
//
//    labelList& m = cellCells.m();
//    const labelList& offsets = cellCells.offsets();
//
//    // 1. For internal faces is just offsetted owner and neighbour
//    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
//    {
//        if (!blockedFace[facei])
//        {
//            label own = agglom[faceOwner[facei]];
//            label nei = agglom[faceNeighbour[facei]];
//
//            m[offsets[own] + nFacesPerCell[own]++] =
//              globalAgglom.toGlobal(nei);
//            m[offsets[nei] + nFacesPerCell[nei]++] =
//              globalAgglom.toGlobal(own);
//        }
//    }
//
//    // 2. For boundary faces is offsetted coupled neighbour
//    forAll(patches, patchi)
//    {
//        const polyPatch& pp = patches[patchi];
//
//        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
//        {
//            label facei = pp.start();
//            label bFacei = pp.start()-mesh.nInternalFaces();
//
//            forAll(pp, i)
//            {
//                if (!blockedFace[facei])
//                {
//                    label own = agglom[faceOwner[facei]];
//
//                    label globalNei = globalNeighbour[bFacei];
//
//                    if
//                    (
//                       !globalAgglom.isLocal(globalNei)
//                     || globalAgglom.toLocal(globalNei) != own
//                    )
//                    {
//                        m[offsets[own] + nFacesPerCell[own]++] = globalNei;
//                    }
//
//                    facei++;
//                    bFacei++;
//                }
//            }
//        }
//    }
//
//    // 3. Explicit connections between non-coupled boundary faces
//    forAll(explicitConnections, i)
//    {
//        const labelPair& baffle = explicitConnections[i];
//        label f0 = baffle.first();
//        label f1 = baffle.second();
//
//        if (!blockedFace[f0] && blockedFace[f1])
//        {
//            label f0Own = agglom[faceOwner[f0]];
//            label f1Own = agglom[faceOwner[f1]];
//
//            // Always count the connection between the two owner sides
//            if (f0Own != f1Own)
//            {
//                m[offsets[f0Own] + nFacesPerCell[f0Own]++] =
//                    globalAgglom.toGlobal(f1Own);
//                m[offsets[f1Own] + nFacesPerCell[f1Own]++] =
//                    globalAgglom.toGlobal(f0Own);
//            }
//
//            // Add any neighbour side connections
//            if (mesh.isInternalFace(f0))
//            {
//                label f0Nei = agglom[faceNeighbour[f0]];
//
//                if (mesh.isInternalFace(f1))
//                {
//                    // Internal faces
//                    label f1Nei = agglom[faceNeighbour[f1]];
//
//                    if (f0Own != f1Nei)
//                    {
//                        m[offsets[f0Own] + nFacesPerCell[f0Own]++] =
//                            globalAgglom.toGlobal(f1Nei);
//                        m[offsets[f1Nei] + nFacesPerCell[f1Nei]++] =
//                            globalAgglom.toGlobal(f1Nei);
//                    }
//                    if (f0Nei != f1Own)
//                    {
//                        m[offsets[f0Nei] + nFacesPerCell[f0Nei]++] =
//                            globalAgglom.toGlobal(f1Own);
//                        m[offsets[f1Own] + nFacesPerCell[f1Own]++] =
//                            globalAgglom.toGlobal(f0Nei);
//                    }
//                    if (f0Nei != f1Nei)
//                    {
//                        m[offsets[f0Nei] + nFacesPerCell[f0Nei]++] =
//                            globalAgglom.toGlobal(f1Nei);
//                        m[offsets[f1Nei] + nFacesPerCell[f1Nei]++] =
//                            globalAgglom.toGlobal(f0Nei);
//                    }
//                }
//                else
//                {
//                    // f1 boundary face
//                    if (f0Nei != f1Own)
//                    {
//                        m[offsets[f0Nei] + nFacesPerCell[f0Nei]++] =
//                            globalAgglom.toGlobal(f1Own);
//                        m[offsets[f1Own] + nFacesPerCell[f1Own]++] =
//                            globalAgglom.toGlobal(f0Nei);
//                    }
//                }
//            }
//            else
//            {
//                if (mesh.isInternalFace(f1))
//                {
//                    label f1Nei = agglom[faceNeighbour[f1]];
//                    if (f0Own != f1Nei)
//                    {
//                        m[offsets[f0Own] + nFacesPerCell[f0Own]++] =
//                            globalAgglom.toGlobal(f1Nei);
//                        m[offsets[f1Nei] + nFacesPerCell[f1Nei]++] =
//                            globalAgglom.toGlobal(f0Own);
//                    }
//                }
//            }
//        }
//    }
//
//
//    // Check for duplicates connections between cells
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    // Done as postprocessing step since we now have cellCells.
//    label newIndex = 0;
//    labelHashSet nbrCells;
//
//
//    if (cellCells.size() == 0)
//    {
//        return;
//    }
//
//    label startIndex = cellCells.offsets()[0];
//
//    forAll(cellCells, celli)
//    {
//        nbrCells.clear();
//        nbrCells.insert(globalAgglom.toGlobal(celli));
//
//        label endIndex = cellCells.offsets()[celli+1];
//
//        for (label i = startIndex; i < endIndex; ++i)
//        {
//            if (nbrCells.insert(cellCells.m()[i]))
//            {
//                cellCells.m()[newIndex++] = cellCells.m()[i];
//            }
//        }
//        startIndex = endIndex;
//        cellCells.offsets()[celli+1] = newIndex;
//    }
//
//    cellCells.m().setSize(newIndex);
//
//    //forAll(cellCells, celli)
//    //{
//    //    Pout<< "Original: Coarse cell " << celli << endl;
//    //    forAll(mesh.cellCells()[celli], i)
//    //    {
//    //        Pout<< "    nbr:" << mesh.cellCells()[celli][i] << endl;
//    //    }
//    //    Pout<< "Compacted: Coarse cell " << celli << endl;
//    //    const labelUList cCells = cellCells[celli];
//    //    forAll(cCells, i)
//    //    {
//    //        Pout<< "    nbr:" << cCells[i] << endl;
//    //    }
//    //}
//}


Foam::labelList Foam::decompositionMethod::decompose
(
    const polyMesh& mesh,
    const scalarField& cellWeights,

    //- Whether owner and neighbour should be on same processor
    //  (takes priority over explicitConnections)
    const boolList& blockedFace,

    //- Whether whole sets of faces (and point neighbours) need to be kept
    //  on single processor
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,

    //- Additional connections between boundary faces
    const List<labelPair>& explicitConnections
)
{
    // Any weights specified?
    label nWeights = returnReduce(cellWeights.size(), sumOp<label>());

    if (nWeights > 0 && cellWeights.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Number of weights " << cellWeights.size()
            << " differs from number of cells " << mesh.nCells()
            << exit(FatalError);
    }


    // Any processor sets?
    label nProcSets = 0;
    forAll(specifiedProcessorFaces, setI)
    {
        nProcSets += specifiedProcessorFaces[setI].size();
    }
    reduce(nProcSets, sumOp<label>());

    // Any non-mesh connections?
    label nConnections = returnReduce
    (
        explicitConnections.size(),
        sumOp<label>()
    );

    // Any faces not blocked?
    label nUnblocked = 0;
    forAll(blockedFace, facei)
    {
        if (!blockedFace[facei])
        {
            nUnblocked++;
        }
    }
    reduce(nUnblocked, sumOp<label>());



    // Either do decomposition on cell centres or on agglomeration

    labelList finalDecomp;


    if (nProcSets+nConnections+nUnblocked == 0)
    {
        // No constraints, possibly weights

        if (nWeights > 0)
        {
            finalDecomp = decompose
            (
                mesh,
                mesh.cellCentres(),
                cellWeights
            );
        }
        else
        {
            finalDecomp = decompose(mesh, mesh.cellCentres());
        }
    }
    else
    {
        if (debug)
        {
            Info<< "Constrained decomposition:" << endl
                << "    faces with same owner and neighbour processor : "
                << nUnblocked << endl
                << "    baffle faces with same owner processor        : "
                << nConnections << endl
                << "    faces all on same processor                   : "
                << nProcSets << endl << endl;
        }

        // Determine local regions, separated by blockedFaces
        regionSplit localRegion(mesh, blockedFace, explicitConnections, false);


        if (debug)
        {
            Info<< "Constrained decomposition:" << endl
                << "    split into " << localRegion.nLocalRegions()
                << " regions."
                << endl;
        }

        // Determine region cell centres
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // This just takes the first cell in the region. Otherwise the problem
        // is with cyclics - if we'd average the region centre might be
        // somewhere in the middle of the domain which might not be anywhere
        // near any of the cells.

        pointField regionCentres(localRegion.nLocalRegions(), point::max);

        forAll(localRegion, celli)
        {
            const label regioni = localRegion[celli];

            if (regionCentres[regioni] == point::max)
            {
                regionCentres[regioni] = mesh.cellCentres()[celli];
            }
        }

        // Do decomposition on agglomeration
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        scalarField regionWeights(localRegion.nLocalRegions(), 0);

        if (nWeights > 0)
        {
            forAll(localRegion, celli)
            {
                const label regioni = localRegion[celli];

                regionWeights[regioni] += cellWeights[celli];
            }
        }
        else
        {
            forAll(localRegion, celli)
            {
                const label regioni = localRegion[celli];

                regionWeights[regioni] += 1.0;
            }
        }

        finalDecomp = decompose
        (
            mesh,
            localRegion,
            regionCentres,
            regionWeights
        );



        // Implement the explicitConnections since above decompose
        // does not know about them
        forAll(explicitConnections, connectioni)
        {
            const labelPair& baffle = explicitConnections[connectioni];
            const label f0 = baffle.first();
            const label f1 = baffle.second();

            if (!blockedFace[f0] && !blockedFace[f1])
            {
                // Note: what if internal faces and owner and neighbour on
                // different processor? So for now just push owner side
                // proc

                const label proci = finalDecomp[mesh.faceOwner()[f0]];

                finalDecomp[mesh.faceOwner()[f1]] = proci;
                if (mesh.isInternalFace(f1))
                {
                    finalDecomp[mesh.faceNeighbour()[f1]] = proci;
                }
            }
            else if (blockedFace[f0] != blockedFace[f1])
            {
                FatalErrorInFunction
                    << "On explicit connection between faces " << f0
                    << " and " << f1
                    << " the two blockedFace status are not equal : "
                    << blockedFace[f0] << " and " << blockedFace[f1]
                    << exit(FatalError);
            }
        }


        // blockedFaces corresponding to processor faces need to be handled
        // separately since not handled by local regionSplit. We need to
        // walk now across coupled faces and make sure to move a whole
        // global region across
        if (Pstream::parRun())
        {
            // Re-do regionSplit

            // Field on cells and faces.
            List<minData> cellData(mesh.nCells());
            List<minData> faceData(mesh.nFaces());

            // Take over blockedFaces by seeding a negative number
            // (so is always less than the decomposition)
            label nUnblocked = 0;
            forAll(blockedFace, facei)
            {
                if (blockedFace[facei])
                {
                    faceData[facei] = minData(-123);
                }
                else
                {
                    nUnblocked++;
                }
            }

            // Seed unblocked faces with destination processor
            labelList seedFaces(nUnblocked);
            List<minData> seedData(nUnblocked);
            nUnblocked = 0;

            forAll(blockedFace, facei)
            {
                if (!blockedFace[facei])
                {
                    const label own = mesh.faceOwner()[facei];
                    seedFaces[nUnblocked] = facei;
                    seedData[nUnblocked] = minData(finalDecomp[own]);
                    nUnblocked++;
                }
            }


            // Propagate information inwards
            FaceCellWave<minData> deltaCalc
            (
                mesh,
                seedFaces,
                seedData,
                faceData,
                cellData,
                mesh.globalData().nTotalCells()+1
            );

            // And extract
            forAll(finalDecomp, celli)
            {
                if (cellData[celli].valid(deltaCalc.data()))
                {
                    finalDecomp[celli] = cellData[celli].data();
                }
            }
        }


        // For specifiedProcessorFaces rework the cellToProc to enforce
        // all on one processor since we can't guarantee that the input
        // to regionSplit was a single region.
        // E.g. faceSet 'a' with the cells split into two regions
        // by a notch formed by two walls
        //
        //          \   /
        //           \ /
        //    ---a----+-----a-----
        //
        //
        // Note that reworking the cellToProc might make the decomposition
        // unbalanced.
        forAll(specifiedProcessorFaces, setI)
        {
            const labelList& set = specifiedProcessorFaces[setI];

            label proci = specifiedProcessor[setI];
            if (proci == -1)
            {
                // If no processor specified use the one from the
                // 0th element
                if (set.size())
                {
                    proci = finalDecomp[mesh.faceOwner()[set[0]]];
                }
                else
                {
                    // Zero-sized processor (e.g. from redistributePar)
                    proci = 0;
                }
            }

            forAll(set, fI)
            {
                const face& f = mesh.faces()[set[fI]];
                forAll(f, fp)
                {
                    const labelList& pFaces = mesh.pointFaces()[f[fp]];
                    for (const label facei : pFaces)
                    {
                        finalDecomp[mesh.faceOwner()[facei]] = proci;
                        if (mesh.isInternalFace(facei))
                        {
                            finalDecomp[mesh.faceNeighbour()[facei]] = proci;
                        }
                    }
                }
            }
        }


        if (debug && Pstream::parRun())
        {
            labelList nbrDecomp;
            syncTools::swapBoundaryCellList(mesh, finalDecomp, nbrDecomp);

            const polyBoundaryMesh& patches = mesh.boundaryMesh();
            forAll(patches, patchi)
            {
                const polyPatch& pp = patches[patchi];
                if (pp.coupled())
                {
                    forAll(pp, i)
                    {
                        label facei = pp.start()+i;
                        label own = mesh.faceOwner()[facei];
                        label bFacei = facei-mesh.nInternalFaces();

                        if (!blockedFace[facei])
                        {
                            label ownProc = finalDecomp[own];
                            label nbrProc = nbrDecomp[bFacei];
                            if (ownProc != nbrProc)
                            {
                                FatalErrorInFunction
                                    << "patch:" << pp.name()
                                    << " face:" << facei
                                    << " at:" << mesh.faceCentres()[facei]
                                    << " ownProc:" << ownProc
                                    << " nbrProc:" << nbrProc
                                    << exit(FatalError);
                            }
                        }
                    }
                }
            }
        }
    }

    return finalDecomp;
}


void Foam::decompositionMethod::setConstraints
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
)
{
    blockedFace.setSize(mesh.nFaces());
    blockedFace = true;

    specifiedProcessorFaces.clear();
    explicitConnections.clear();

    forAll(constraints_, constrainti)
    {
        constraints_[constrainti].add
        (
            mesh,
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            explicitConnections
        );
    }
}


void Foam::decompositionMethod::applyConstraints
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
)
{
    forAll(constraints_, constrainti)
    {
        constraints_[constrainti].apply
        (
            mesh,
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            explicitConnections,
            decomposition
        );
    }
}


Foam::labelList Foam::decompositionMethod::decompose
(
    const polyMesh& mesh,
    const scalarField& cellWeights
)
{
    // Collect all constraints

    boolList blockedFace;
    PtrList<labelList> specifiedProcessorFaces;
    labelList specifiedProcessor;
    List<labelPair> explicitConnections;
    setConstraints
    (
        mesh,
        blockedFace,
        specifiedProcessorFaces,
        specifiedProcessor,
        explicitConnections
    );


    // Construct decomposition method and either do decomposition on
    // cell centres or on agglomeration

    labelList finalDecomp = decompose
    (
        mesh,
        cellWeights,            // optional weights
        blockedFace,            // any cells to be combined
        specifiedProcessorFaces,// any whole cluster of cells to be kept
        specifiedProcessor,
        explicitConnections     // baffles
    );


    // Give any constraint the option of modifying the decomposition

    applyConstraints
    (
        mesh,
        blockedFace,
        specifiedProcessorFaces,
        specifiedProcessor,
        explicitConnections,
        finalDecomp
    );

    return finalDecomp;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "multiLevelDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "globalIndex.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiLevelDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        multiLevelDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Given a subset of cells determine the new global indices. The problem
// is in the cells from neighbouring processors which need to be renumbered.
void Foam::multiLevelDecomp::subsetGlobalCellCells
(
    const label nDomains,
    const label domainI,
    const labelList& dist,

    const labelListList& cellCells,
    const labelList& set,
    labelListList& subCellCells,
    labelList& cutConnections
) const
{
    // Determine new index for cells by inverting subset
    labelList oldToNew(invert(cellCells.size(), set));

    globalIndex globalCells(cellCells.size());

    // Subset locally the elements for which I have data
    subCellCells = UIndirectList<labelList>(cellCells, set);

    // Get new indices for neighbouring processors
    List<Map<label>> compactMap;
    mapDistribute map(globalCells, subCellCells, compactMap);
    map.distribute(oldToNew);
    labelList allDist(dist);
    map.distribute(allDist);

    // Now we have:
    // oldToNew : the locally-compact numbering of all our cellCells. -1 if
    //            cellCell is not in set.
    // allDist  : destination domain for all our cellCells
    // subCellCells : indexes into oldToNew and allDist

    // Globally compact numbering for cells in set.
    globalIndex globalSubCells(set.size());

    // Now subCellCells contains indices into oldToNew which are the
    // new locations of the neighbouring cells.

    cutConnections.setSize(nDomains);
    cutConnections = 0;

    forAll(subCellCells, subCelli)
    {
        labelList& cCells = subCellCells[subCelli];

        // Keep the connections to valid mapped cells
        label newI = 0;
        forAll(cCells, i)
        {
            // Get locally-compact cell index of neighbouring cell
            const label nbrCelli = oldToNew[cCells[i]];
            if (nbrCelli == -1)
            {
                cutConnections[allDist[cCells[i]]]++;
            }
            else
            {
                // Reconvert local cell index into global one

                // Get original neighbour
                const label celli = set[subCelli];
                const label oldNbrCelli = cellCells[celli][i];
                // Get processor from original neighbour
                const label proci = globalCells.whichProcID(oldNbrCelli);
                // Convert into global compact numbering
                cCells[newI++] = globalSubCells.toGlobal(proci, nbrCelli);
            }
        }
        cCells.setSize(newI);
    }
}


void Foam::multiLevelDecomp::decompose
(
    const labelListList& pointPoints,
    const pointField& points,
    const scalarField& pointWeights,
    const labelUList& pointMap,     // map back to original points
    const label currLevel,
    const label leafOffset,

    labelList& finalDecomp
)
{
    labelList dist
    (
        methods_[currLevel].decompose
        (
            pointPoints,
            points,
            pointWeights
        )
    );

    // The next recursion level
    const label nextLevel = currLevel+1;

    // Number of domains at this current level
    const label nCurrDomains = methods_[currLevel].nDomains();

    // Calculate the domain remapping.
    // The decompose() method delivers a distribution of [0..nDomains-1]
    // which we map to the final location according to the decomposition
    // leaf we are on.

    labelList domainLookup(nCurrDomains);
    {
        label sizes = 1;  // Cumulative number of domains
        for (label i = 0; i <= currLevel; ++i)
        {
            sizes *= methods_[i].nDomains();
        }

        // Distribution of domains at this level
        sizes = this->nDomains() / sizes;

        forAll(domainLookup, i)
        {
            domainLookup[i] = i * sizes + leafOffset;
        }
    }

    if (debug)
    {
        Info<< "Distribute at level " << currLevel
            << " to domains" << nl
            << flatOutput(domainLookup) << endl;
    }

    // Extract processor+local index from point-point addressing
    forAll(pointMap, i)
    {
        const label orig = pointMap[i];
        finalDecomp[orig] = domainLookup[dist[i]];
    }

    if (nextLevel < methods_.size())
    {
        // Recurse

        // Determine points per domain
        labelListList domainToPoints(invertOneToMany(nCurrDomains, dist));

        // Extract processor+local index from point-point addressing
        if (debug && Pstream::master())
        {
            Pout<< "Decomposition at level " << currLevel << " :" << endl;
        }

        for (label domainI = 0; domainI < nCurrDomains; domainI++)
        {
            // Extract elements for current domain
            const labelList domainPoints(findIndices(dist, domainI));

            // Subset point-wise data.
            pointField subPoints(points, domainPoints);
            scalarField subWeights(pointWeights, domainPoints);
            labelList subPointMap(UIndirectList<label>(pointMap, domainPoints));
            // Subset point-point addressing (adapt global numbering)
            labelListList subPointPoints;
            labelList nOutsideConnections;
            subsetGlobalCellCells
            (
                nCurrDomains,
                domainI,
                dist,

                pointPoints,
                domainPoints,

                subPointPoints,
                nOutsideConnections
            );

            label nPoints = returnReduce(domainPoints.size(), plusOp<label>());
            Pstream::listCombineGather(nOutsideConnections, plusEqOp<label>());
            Pstream::listCombineScatter(nOutsideConnections);
            label nPatches = 0;
            label nFaces = 0;
            for (const label nConnect : nOutsideConnections)
            {
                if (nConnect > 0)
                {
                    ++nPatches;
                    nFaces += nConnect;
                }
            }

            string oldPrefix;
            if (debug && Pstream::master())
            {
                Pout<< "    Domain " << domainI << nl
                    << "        Number of cells = " << nPoints << nl
                    << "        Number of inter-domain patches = " << nPatches
                    << nl
                    << "        Number of inter-domain faces = " << nFaces << nl
                    << endl;
                oldPrefix = Pout.prefix();
                Pout.prefix() = "  " + oldPrefix;
            }

            decompose
            (
                subPointPoints,
                subPoints,
                subWeights,
                subPointMap,
                nextLevel,
                domainLookup[domainI], // The offset for this level and leaf

                finalDecomp
            );
            if (debug && Pstream::master())
            {
                Pout.prefix() = oldPrefix;
            }
        }


        if (debug)
        {
            // Do straight decompose of two levels
            const label nNext = methods_[nextLevel].nDomains();
            const label nTotal = nCurrDomains * nNext;

            // Retrieve original level0 dictionary and modify number of domains
            dictionary::const_iterator iter =
                decompositionDict_.optionalSubDict(typeName + "Coeffs").begin();
            dictionary myDict = iter().dict();
            myDict.set("numberOfSubdomains", nTotal);

            if (debug && Pstream::master())
            {
                Pout<< "Reference decomposition with " << myDict << " :"
                    << endl;
            }

            autoPtr<decompositionMethod> method0 = decompositionMethod::New
            (
                myDict
            );
            labelList dist
            (
                method0().decompose
                (
                    pointPoints,
                    points,
                    pointWeights
                )
            );

            for (label blockI = 0; blockI < nCurrDomains; blockI++)
            {
                // Count the number inbetween blocks of nNext size

                label nPoints = 0;
                labelList nOutsideConnections(nCurrDomains, 0);
                forAll(pointPoints, pointi)
                {
                    if ((dist[pointi] / nNext) == blockI)
                    {
                        nPoints++;

                        const labelList& pPoints = pointPoints[pointi];

                        forAll(pPoints, i)
                        {
                            const label distBlockI = dist[pPoints[i]] / nNext;
                            if (distBlockI != blockI)
                            {
                                nOutsideConnections[distBlockI]++;
                            }
                        }
                    }
                }

                reduce(nPoints, plusOp<label>());
                Pstream::listCombineGather
                (
                    nOutsideConnections,
                    plusEqOp<label>()
                );
                Pstream::listCombineScatter(nOutsideConnections);
                label nPatches = 0;
                label nFaces = 0;
                for (const label nConnect : nOutsideConnections)
                {
                    if (nConnect > 0)
                    {
                        ++nPatches;
                        nFaces += nConnect;
                    }
                }

                if (debug && Pstream::master())
                {
                    Pout<< "    Domain " << blockI << nl
                        << "        Number of cells = " << nPoints << nl
                        << "        Number of inter-domain patches = "
                        << nPatches << nl
                        << "        Number of inter-domain faces = " << nFaces
                        << nl << endl;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiLevelDecomp::multiLevelDecomp(const dictionary& decompositionDict)
:
    decompositionMethod(decompositionDict),
    methodsDict_(decompositionDict_.optionalSubDict(typeName + "Coeffs"))
{
    methods_.setSize(methodsDict_.size());
    label i = 0;
    forAllConstIter(dictionary, methodsDict_, iter)
    {
        methods_.set(i++, decompositionMethod::New(iter().dict()));
    }

    label n = 1;
    Info<< "decompositionMethod " << type() << " :" << endl;
    forAll(methods_, i)
    {
        Info<< "    level " << i << " decomposing with " << methods_[i].type()
            << " into " << methods_[i].nDomains() << " subdomains." << endl;

        n *= methods_[i].nDomains();
    }

    if (n != nDomains())
    {
        FatalErrorInFunction
            << "Top level decomposition specifies " << nDomains()
            << " domains which is not equal to the product of"
            << " all sub domains " << n
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiLevelDecomp::parallelAware() const
{
    forAll(methods_, i)
    {
        if (!methods_[i].parallelAware())
        {
            return false;
        }
    }
    return true;
}


Foam::labelList Foam::multiLevelDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& cc,
    const scalarField& cWeights
)
{
    CompactListList<label> cellCells;
    calcCellCells(mesh, identity(cc.size()), cc.size(), true, cellCells);

    labelList finalDecomp(cc.size(), 0);
    labelList cellMap(identity(cc.size()));

    decompose
    (
        cellCells(),
        cc,
        cWeights,
        cellMap,      // map back to original cells
        0,
        0,

        finalDecomp
    );

    return finalDecomp;
}


Foam::labelList Foam::multiLevelDecomp::decompose
(
    const labelListList& globalPointPoints,
    const pointField& points,
    const scalarField& pointWeights
)
{
    labelList finalDecomp(points.size(), 0);
    labelList pointMap(identity(points.size()));

    decompose
    (
        globalPointPoints,
        points,
        pointWeights,
        pointMap,       // map back to original points
        0,
        0,

        finalDecomp
    );

    return finalDecomp;
}


// ************************************************************************* //

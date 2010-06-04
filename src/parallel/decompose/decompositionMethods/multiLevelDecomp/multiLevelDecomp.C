/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    const labelListList& cellCells,
    const labelList& set,
    labelListList& subCellCells,
    label& cutConnections
) const
{
    // Determine new index for cells by inverting subset
    labelList oldToNew(invert(cellCells.size(), set));

    globalIndex globalCells(cellCells.size());

    // Subset locally the elements for which I have data
    subCellCells = UIndirectList<labelList>(cellCells, set);

    // Get new indices for neighbouring processors
    List<Map<label> > compactMap;
    mapDistribute map(globalCells, subCellCells, compactMap);
    map.distribute(oldToNew);

    // Now subCellCells contains indices into oldToNew which are the
    // new locations of the neighbouring cells.

    cutConnections = 0;

    forAll(subCellCells, subCellI)
    {
        labelList& cCells = subCellCells[subCellI];

        // Keep the connections to valid mapped cells
        label newI = 0;
        forAll(cCells, i)
        {
            label subCellI = oldToNew[cCells[i]];
            if (subCellI == -1)
            {
                cutConnections++;
            }
            else
            {
                cCells[newI++] = subCellI;
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
    const labelList& pointMap,      // map back to original points
    const label levelI,

    labelField& finalDecomp
)
{
    labelList dist
    (
        methods_[levelI].decompose
        (
            pointPoints,
            points,
            pointWeights
        )
    );

//Pout<< "At level " << levelI << endl;
//forAll(dist, i)
//{
//    Pout<< "    " << i << " at:" << points[i]
//        << " original:" << pointMap[i] << "  dist:" << dist[i]
//        << " connected to:" << pointField(points, pointPoints[i])
//        << endl;
//}
//Pout<< endl;

    forAll(pointMap, i)
    {
        label orig = pointMap[i];
        finalDecomp[orig] += dist[i];
    }

    if (levelI != methods_.size()-1)
    {
        // Recurse

        // Determine points per domain
        label n = methods_[levelI].nDomains();
        labelListList domainToPoints(invertOneToMany(n, dist));

        // 'Make space' for new levels of decomposition
        finalDecomp *= methods_[levelI+1].nDomains();

        // Extract processor+local index from point-point addressing

        forAll(domainToPoints, domainI)
        {
            const labelList& myPoints = domainToPoints[domainI];

            // Subset point-wise data.
            pointField subPoints(points, myPoints);
            scalarField subWeights(pointWeights, myPoints);
            labelList subPointMap(UIndirectList<label>(pointMap, myPoints));
            // Subset point-point addressing (adapt global numbering)
            labelListList subPointPoints;
            label nOutsideConnections;
            subsetGlobalCellCells
            (
                pointPoints,
                myPoints,
                subPointPoints,
                nOutsideConnections
            );

            //Info<< "At level " << levelI << "  domain " << domainI
            //    << " have connections to other domains "
            //    << returnReduce(nOutsideConnections, sumOp<label>())
            //    << endl;

            decompose
            (
                subPointPoints,
                subPoints,
                subWeights,
                subPointMap,
                levelI+1,

                finalDecomp
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiLevelDecomp::multiLevelDecomp(const dictionary& decompositionDict)
:
    decompositionMethod(decompositionDict)
{
    const dictionary& myDict = decompositionDict_.subDict(typeName + "Coeffs");

    methods_.setSize(myDict.size());
    label i = 0;
    forAllConstIter(dictionary, myDict, iter)
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
        FatalErrorIn("multiLevelDecomp::multiLevelDecomp(const dictionary&)")
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
    calcCellCells(mesh, identity(cc.size()), cc.size(), cellCells);

    labelField finalDecomp(cc.size(), 0);
    labelList cellMap(identity(cc.size()));

    decompose
    (
        cellCells(),
        cc,
        cWeights,
        cellMap,      // map back to original cells
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
    labelField finalDecomp(points.size(), 0);
    labelList pointMap(identity(points.size()));

    decompose
    (
        globalPointPoints,
        points,
        pointWeights,
        pointMap,       // map back to original points
        0,

        finalDecomp
    );

    return finalDecomp;
}


// ************************************************************************* //

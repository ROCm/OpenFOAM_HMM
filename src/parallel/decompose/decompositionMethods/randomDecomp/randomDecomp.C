/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "randomDecomp.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(randomDecomp);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        randomDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::randomDecomp::randomMap(const label nCells) const
{
    Random rndGen(0);

    labelList finalDecomp(nCells);

    for (label& val : finalDecomp)
    {
        val = rndGen.position<label>(0, nDomains_ - 1);
    }

    return finalDecomp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::randomDecomp::randomDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    decompositionMethod(decompDict, regionName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::randomDecomp::decompose
(
    const pointField& points,
    const scalarField& pointWeights
) const
{
    return randomMap(points.size());
}


Foam::labelList Foam::randomDecomp::decompose
(
    const pointField& points
) const
{
    return randomMap(points.size());
}


Foam::labelList Foam::randomDecomp::decompose
(
    const polyMesh& mesh,
    const pointField&,
    const scalarField&
) const
{
    return randomMap(mesh.nCells());  // or cc.size()
}


Foam::labelList Foam::randomDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField&,
    const scalarField&
) const
{
    return randomMap(globalCellCells.size());  // or cc.size()
}


// ************************************************************************* //

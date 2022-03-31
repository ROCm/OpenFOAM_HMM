/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "springRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "decompositionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(springRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        springRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::springRenumber::springRenumber(const dictionary& dict)
:
    renumberMethod(dict),
    coeffsDict_(dict.optionalSubDict(typeName+"Coeffs")),
    maxIter_(coeffsDict_.get<label>("maxIter")),
    maxCo_(coeffsDict_.get<scalar>("maxCo")),
    freezeFraction_(coeffsDict_.get<scalar>("freezeFraction"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ConnectionListListType>
Foam::labelList Foam::springRenumber::renumberImpl
(
    const ConnectionListListType& cellCells
) const
{
    const label nOldCells(cellCells.size());

    // Look at cell index as a 1D position parameter.
    // Move cells to the average 'position' of their neighbour.

    scalarField position(nOldCells);
    forAll(position, celli)
    {
        position[celli] = celli;
    }

    // Sum force per cell. Also reused for the displacement
    scalarField sumForce(nOldCells);

    labelList oldToNew(identity(nOldCells));

    scalar maxCo = (maxCo_ * nOldCells);

    for (label iter = 0; iter < maxIter_; ++iter)
    {
        //Pout<< "Iteration : " << iter << nl
        //    << "------------"
        //    << endl;

        //Pout<< "Position :" << nl
        //    << "    min : " << min(position) << nl
        //    << "    max : " << max(position) << nl
        //    << "    avg : " << average(position) << nl
        //    << endl;

        // Sum force per cell.
        sumForce = Zero;
        for (label oldCelli = 0; oldCelli < nOldCells; ++oldCelli)
        {
            const label celli = oldToNew[oldCelli];
            const auto& neighbours = cellCells[oldCelli];

            for (const label nbr : neighbours)
            {
                const label nbrCelli = oldToNew[nbr];

                sumForce[celli] += (position[nbrCelli]-position[celli]);
            }
        }

        //Pout<< "Force :" << nl
        //    << "    min    : " << min(sumForce) << nl
        //    << "    max    : " << max(sumForce) << nl
        //    << "    avgMag : " << average(mag(sumForce)) << nl
        //    << "DeltaT : " << deltaT << nl
        //    << endl;

        // Limit displacement
        scalar deltaT = maxCo / max(mag(sumForce));

        Info<< "Iter:" << iter
            << "  maxCo:" << maxCo
            << "  deltaT:" << deltaT
            << "  average force:" << average(mag(sumForce)) << endl;

        // Determine displacement
        scalarField& displacement = sumForce;
        displacement *= deltaT;

        //Pout<< "Displacement :" << nl
        //    << "    min    : " << min(displacement) << nl
        //    << "    max    : " << max(displacement) << nl
        //    << "    avgMag : " << average(mag(displacement)) << nl
        //    << endl;

        // Calculate new position and scale to be within original range
        // (0..nCells-1) for ease of postprocessing.
        position += displacement;
        position -= min(position);
        position *= (position.size()-1)/max(position);

        // Slowly freeze.
        maxCo *= freezeFraction_;
    }

    //writeOBJ("endPosition.obj", cellCells, position);

    // Move cells to new position
    labelList shuffle(sortedOrder(position));

    // Reorder oldToNew
    inplaceReorder(shuffle, oldToNew);

    return invert(nOldCells, oldToNew);
}


Foam::labelList Foam::springRenumber::renumber
(
    const polyMesh& mesh,
    const pointField&
) const
{
    CompactListList<label> cellCells;
    decompositionMethod::calcCellCells
    (
        mesh,
        identity(mesh.nCells()),
        mesh.nCells(),
        false,                      // local only
        cellCells
    );

    return renumberImpl(cellCells);
}


Foam::labelList Foam::springRenumber::renumber
(
    const CompactListList<label>& cellCells,
    const pointField&
) const
{
    return renumberImpl(cellCells);
}


Foam::labelList Foam::springRenumber::renumber
(
    const labelListList& cellCells,
    const pointField&
) const
{
    return renumberImpl(cellCells);
}


// ************************************************************************* //

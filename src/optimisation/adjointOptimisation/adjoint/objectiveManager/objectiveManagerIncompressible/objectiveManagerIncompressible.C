/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "objectiveManagerIncompressible.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveManagerIncompressible, 0);
addToRunTimeSelectionTable
(
    objectiveManager,
    objectiveManagerIncompressible,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveManagerIncompressible::objectiveManagerIncompressible
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveManager(mesh, dict, adjointSolverName, primalSolverName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void objectiveManagerIncompressible::addUaEqnSource(fvVectorMatrix& UaEqn)
{
    // Add contributions from objective functions
    for (objective& obj : objectives_)
    {
        auto& icoObj = refCast<objectiveIncompressible>(obj);

        if (icoObj.hasdJdv())
        {
            scalar weight = icoObj.weight();
            UaEqn += weight*icoObj.dJdv();
        }
    }
}


void objectiveManagerIncompressible::addPaEqnSource(fvScalarMatrix& paEqn)
{
    // Add contributions from objective functions
    for (objective& obj : objectives_)
    {
        auto& icoObj = refCast<objectiveIncompressible>(obj);

        if (icoObj.hasdJdp())
        {
            scalar weight = icoObj.weight();
            paEqn += weight*icoObj.dJdp();
        }
    }
}


void objectiveManagerIncompressible::addTMEqn1Source(fvScalarMatrix& adjTMEqn1)
{
    // Add contributions from objective functions
    for (objective& obj : objectives_)
    {
        auto& icoObj = refCast<objectiveIncompressible>(obj);

        if (icoObj.hasdJdTMVar1())
        {
            scalar weight = icoObj.weight();
            adjTMEqn1 += weight*icoObj.dJdTMvar1();
        }
    }
}


void objectiveManagerIncompressible::addTMEqn2Source(fvScalarMatrix& adjTMEqn2)
{
    // Add contributions from objective functions
    for (objective& obj : objectives_)
    {
        auto& icoObj = refCast<objectiveIncompressible>(obj);

        if (icoObj.hasdJdTMVar2())
        {
            scalar weight = icoObj.weight();
            adjTMEqn2 += weight*icoObj.dJdTMvar2();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

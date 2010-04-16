/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "fvData.H"
#include "Time.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::fvData::debug(Foam::debug::debugSwitch("fvData", false));

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvData::fvData(const objectRegistry& obr)
:
    IOdictionary
    (
        IOobject
        (
            "fvData",
            obr.time().system(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    )
{
    set("solverPerformance", dictionary());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::fvData::solverPerformanceDict() const
{
    return subDict("solverPerformance");
}


void Foam::fvData::setSolverPerformance
(
    const word& name,
    const lduMatrix::solverPerformance& sp
) const
{
    const_cast<dictionary&>(solverPerformanceDict()).set(name, sp);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "solution.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int Foam::solution::debug(Foam::debug::debugSwitch("solution", false));

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solution::solution(const objectRegistry& obr, const fileName& dictName)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr.time().system(),
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    relaxationFactors_(ITstream("relaxationFactors", tokenList())()),
    solvers_(ITstream("solvers", tokenList())())
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solution::read()
{
    if (regIOobject::read())
    {
        const dictionary& dict = solutionDict();

        if (dict.found("relaxationFactors"))
        {
            relaxationFactors_ = dict.subDict("relaxationFactors");
        }

        if (dict.found("solvers"))
        {
            solvers_ = dict.subDict("solvers");
        }

        return true;
    }
    else
    {
        return false;
    }
}


const Foam::dictionary& Foam::solution::solutionDict() const
{
    if (found("select"))
    {
        return subDict(word(lookup("select")));
    }
    else
    {
        return *this;
    }
}


bool Foam::solution::relax(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup relax for " << name << endl;
    }

    return relaxationFactors_.found(name);
}


Foam::scalar Foam::solution::relaxationFactor(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup relaxationFactor for " << name << endl;
    }

    return readScalar(relaxationFactors_.lookup(name));
}


const Foam::dictionary& Foam::solution::solverDict(const word& name) const
{
    if (debug)
    {
        InfoIn("solution::solverDict(const word& name)")
            << "Lookup solver for " << name << endl;
    }

    return solvers_.subDict(name);
}


Foam::ITstream& Foam::solution::solver(const word& name) const
{
    if (debug)
    {
        InfoIn("solution::solver(const word& name)")
            << "Lookup solver for " << name << endl;
    }

    return solvers_.lookup(name);
}


// ************************************************************************* //

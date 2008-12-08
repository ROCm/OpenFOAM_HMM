/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

// these two are for old syntax compatibility:
#include "BICCG.H"
#include "ICCG.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int Foam::solution::debug(::Foam::debug::debugSwitch("solution", 0));

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

Foam::label Foam::solution::upgradeSolverDict
(
    dictionary& dict,
    const bool verbose
)
{
    label nChanged = 0;

    // backward compatibility:
    // recast primitive entries into dictionary entries
    forAllIter(dictionary, dict, iter)
    {
        if (!iter().isDict())
        {
            Istream& is = iter().stream();
            word name(is);
            dictionary subdict;

            if (name == "BICCG")
            {
                // special treatment for very old syntax
                subdict = BICCG::solverDict(is);
            }
            else if (name == "ICCG")
            {
                // special treatment for very old syntax
                subdict = ICCG::solverDict(is);
            }
            else
            {
                subdict.add("solver", name);
                subdict <<= dictionary(is);

                // preconditioner can be a primitiveEntry w/o settings,
                // or a dictionaryEntry.
                // transform primitiveEntry with settings -> dictionaryEntry
                entry* precond = subdict.lookupEntryPtr
                (
                    "preconditioner",
                    false,
                    false
                );

                if (precond && !precond->isDict())
                {
                    Istream& is = precond->stream();
                    is >> name;

                    if (!is.eof())
                    {
                        dictionary precondDict;
                        precondDict.add("preconditioner", name);
                        precondDict <<= dictionary(is);

                        subdict.set("preconditioner", precondDict);
                    }
                }
            }

            // write out information to help people adjust to the new syntax
            if (verbose)
            {
                Info<< "// using new solver syntax:\n"
                    << iter().keyword() << subdict << endl;
            }

            // overwrite with dictionary entry
            dict.set(iter().keyword(), subdict);

            nChanged++;
        }
    }

    return nChanged;
}


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
            upgradeSolverDict(solvers_);
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
        Info<< "Find relax for " << name << endl;
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


const Foam::dictionary& Foam::solution::solver(const word& name) const
{
    if (debug)
    {
        InfoIn("solution::solver(const word& name)")
            << "Lookup solver for " << name << endl;
    }

    return solvers_.subDict(name);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "solution.H"
#include "HashPtrTable.H"
#include "Function1.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineDebugSwitchWithName(solution, "solution", 0);
}

// List of sub-dictionaries to rewrite
static const Foam::List<Foam::word> subDictNames
({
    "preconditioner", "smoother"
});


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::solution::read(const dictionary& dict)
{
    if (dict.found("cache"))
    {
        cache_ = dict.subDict("cache");
        caching_ = cache_.getOrDefault("active", true);
    }

    if (dict.found("relaxationFactors"))
    {
        const dictionary& relaxDict = dict.subDict("relaxationFactors");

        if (relaxDict.found("fields") || relaxDict.found("equations"))
        {
            if (relaxDict.found("fields"))
            {
                fieldRelaxDict_ = relaxDict.subDict("fields");
                fieldRelaxCache_.clear();
            }

            if (relaxDict.found("equations"))
            {
                eqnRelaxDict_ = relaxDict.subDict("equations");
                eqnRelaxCache_.clear();
            }
        }
        else
        {
            // backwards compatibility
            fieldRelaxDict_.clear();
            fieldRelaxCache_.clear();

            for (const word& e : relaxDict.toc())
            {
                scalar value = relaxDict.get<scalar>(e);

                if (e.starts_with('p'))
                {
                    fieldRelaxDict_.add(e, value);
                }
                else if (e.starts_with("rho"))
                {
                    fieldRelaxDict_.add(e, value);
                }
            }

            eqnRelaxDict_ = relaxDict;
            eqnRelaxCache_.clear();
        }


        fieldRelaxDefault_ = Function1<scalar>::NewIfPresent
        (
            "default",
            fieldRelaxDict_
        );
        if (!fieldRelaxDefault_)
        {
            fieldRelaxDefault_.reset
            (
                new Function1Types::Constant<scalar>("default", 0)
            );
        }

        eqnRelaxDefault_ = Function1<scalar>::NewIfPresent
        (
            "default",
            eqnRelaxDict_
        );
        if (!eqnRelaxDefault_)
        {
            eqnRelaxDefault_.reset
            (
                new Function1Types::Constant<scalar>("default", 0)
            );
        }

        DebugInfo
            << "Relaxation factors:" << nl
            << "fields: " << fieldRelaxDict_ << nl
            << "equations: " << eqnRelaxDict_ << endl;
    }

    if (dict.found("solvers"))
    {
        solvers_ = dict.subDict("solvers");
        upgradeSolverDict(solvers_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solution::solution
(
    const objectRegistry& obr,
    const fileName& dictName,
    const dictionary* fallback
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr.time().system(),
            obr,
            (
                obr.readOpt() == IOobject::MUST_READ
             || obr.readOpt() == IOobject::READ_IF_PRESENT
              ? IOobject::MUST_READ_IF_MODIFIED
              : obr.readOpt()
            ),
            IOobject::NO_WRITE
        ),
        fallback
    ),
    cache_(),
    caching_(false),
    fieldRelaxDict_(),
    eqnRelaxDict_(),
    solvers_()
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read(solutionDict());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// A non-default destructor since we had incomplete types in the header
Foam::solution::~solution()
{}


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
    for (const entry& dEntry : dict)
    {
        if (!dEntry.isDict())
        {
            ITstream& is = dEntry.stream();
            word name(is);
            dictionary subdict;

            subdict.add("solver", name);
            subdict <<= dictionary(is);

            // preconditioner and smoother entries can be
            // 1) primitiveEntry w/o settings,
            // 2) or a dictionaryEntry.
            // transform primitiveEntry with settings -> dictionaryEntry
            for (const word& dictName : subDictNames)
            {
                entry* eptr = subdict.findEntry(dictName, keyType::LITERAL);

                if (eptr && !eptr->isDict())
                {
                    ITstream& is = eptr->stream();
                    is >> name;

                    if (!is.eof())
                    {
                        dictionary newDict;
                        newDict.add(dictName, name);
                        newDict <<= dictionary(is);

                        subdict.set(dictName, newDict);
                    }
                }
            }

            // write out information to help people adjust to the new syntax
            if (verbose && Pstream::master())
            {
                Info<< "// using new solver syntax:\n"
                    << dEntry.keyword() << subdict << endl;
            }

            // overwrite with dictionary entry
            dict.set(dEntry.keyword(), subdict);

            ++nChanged;
        }
    }

    return nChanged;
}


bool Foam::solution::cache(const word& name) const
{
    if (caching_)
    {
        DebugInfo<< "Cache: find entry for " << name << endl;
        return cache_.found(name);
    }

    return false;
}


bool Foam::solution::relaxField(const word& name) const
{
    DebugInfo
        << "Field relaxation factor for " << name
        << " is " << (fieldRelaxDict_.found(name) ? "set" : "unset") << endl;

    return fieldRelaxDict_.found(name) || fieldRelaxDict_.found("default");
}


bool Foam::solution::relaxEquation(const word& name) const
{
    DebugInfo<< "Find equation relaxation factor for " << name << endl;
    return eqnRelaxDict_.found(name) || eqnRelaxDict_.found("default");
}


Foam::scalar Foam::solution::fieldRelaxationFactor(const word& name) const
{
    DebugInfo<< "Lookup variable relaxation factor for " << name << endl;

    if (fieldRelaxDict_.found(name))
    {
        return Function1<scalar>::New
        (
            fieldRelaxCache_,  // cache
            name,
            fieldRelaxDict_,
            keyType::REGEX
        )().value(time().timeOutputValue());
    }
    else if (fieldRelaxDefault_)
    {
        return fieldRelaxDefault_().value(time().timeOutputValue());
    }

    FatalIOErrorInFunction(fieldRelaxDict_)
        << "Cannot find variable relaxation factor for '" << name
        << "' or a suitable default value." << nl
        << exit(FatalIOError);

    return 0;
}


Foam::scalar Foam::solution::equationRelaxationFactor(const word& name) const
{
    DebugInfo<< "Lookup equation relaxation factor for " << name << endl;

    if (eqnRelaxDict_.found(name))
    {
        return Function1<scalar>::New
        (
            eqnRelaxCache_,  // cache
            name,
            eqnRelaxDict_,
            keyType::REGEX
        )().value(time().timeOutputValue());
    }
    else if (eqnRelaxDefault_)
    {
        return eqnRelaxDefault_().value(time().timeOutputValue());
    }

    FatalIOErrorInFunction(eqnRelaxDict_)
        << "Cannot find equation relaxation factor for '" << name
        << "' or a suitable default value."
        << exit(FatalIOError);

    return 0;
}


const Foam::dictionary& Foam::solution::solutionDict() const
{
    if (found("select"))
    {
        return subDict(get<word>("select"));
    }

    return *this;
}


const Foam::dictionary& Foam::solution::solverDict(const word& name) const
{
    DebugInfo<< "Lookup solver for " << name << endl;
    return solvers_.subDict(name);
}


const Foam::dictionary& Foam::solution::solver(const word& name) const
{
    DebugInfo<< "Lookup solver for " << name << endl;
    return solvers_.subDict(name);
}


bool Foam::solution::read()
{
    if (regIOobject::read())
    {
        read(solutionDict());

        return true;
    }

    return false;
}


// ************************************************************************* //

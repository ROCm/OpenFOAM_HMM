/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "topoSetSource.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "bitSet.H"
#include "topoSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topoSetSource, 0);
    defineRunTimeSelectionTable(topoSetSource, word);
    defineRunTimeSelectionTable(topoSetSource, istream);
}


Foam::HashTable<Foam::string>* Foam::topoSetSource::usageTablePtr_ = nullptr;


const Foam::Enum
<
    Foam::topoSetSource::setAction
>
Foam::topoSetSource::actionNames
({
    { setAction::ADD, "add" },
    { setAction::SUBTRACT, "subtract" },
    { setAction::SUBSET, "subset" },
    { setAction::INVERT, "invert" },
    { setAction::CLEAR, "clear" },
    { setAction::NEW, "new" },
    { setAction::REMOVE, "remove" },
    { setAction::LIST, "list" },
    { setAction::SUBTRACT, "delete" },   // Compat (1806)
});


const Foam::string Foam::topoSetSource::illegalSource_
(
    "Illegal topoSetSource name"
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::topoSetSource::check(labelList& list, const label maxLabel)
{
    const label len = list.size();

    label nGood = 0;

    for (label i=0; i < len; ++i)
    {
        const label val = list[i];

        if (val >= 0 && val < maxLabel)
        {
            if (nGood != i)
            {
                list[nGood] = val;
            }
            ++nGood;
        }
    }

    const label nReject = (len - nGood);

    if (nReject)
    {
        list.resize(nGood);

        // Report?
    }

    return !nReject;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::topoSetSource> Foam::topoSetSource::New
(
    const word& topoSetSourceType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    auto* ctorPtr = wordConstructorTable(topoSetSourceType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "topoSetSource",
            topoSetSourceType,
            *wordConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<topoSetSource>(ctorPtr(mesh, dict));
}


Foam::autoPtr<Foam::topoSetSource> Foam::topoSetSource::New
(
    const word& topoSetSourceType,
    const polyMesh& mesh,
    Istream& is
)
{
    auto* ctorPtr = istreamConstructorTable(topoSetSourceType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "topoSetSource",
            topoSetSourceType,
            *istreamConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<topoSetSource>(ctorPtr(mesh, is));
}


Foam::Istream& Foam::topoSetSource::checkIs(Istream& is)
{
    if (!is.good() || is.eof())
    {
        FatalErrorInFunction
            << exit(FatalError);
    }

    return is;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::topoSetSource::addOrDelete
(
    topoSet& set,
    const label id,
    const bool add
) const
{
    if (add)
    {
        set.set(id);
    }
    else
    {
        set.unset(id);
    }
}


void Foam::topoSetSource::addOrDelete
(
    topoSet& set,
    const labelUList& labels,
    const bool add
) const
{
    if (add)
    {
        set.set(labels);
    }
    else
    {
        set.unset(labels);
    }
}


void Foam::topoSetSource::addOrDelete
(
    topoSet& set,
    const bitSet& labels,
    const bool add
) const
{
    if (add)
    {
        for (const label id : labels)
        {
            set.set(id);
        }
    }
    else
    {
        for (const label id : labels)
        {
            set.unset(id);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoSetSource::topoSetSource
(
    const polyMesh& mesh,
    bool verbose
)
:
    mesh_(mesh),
    verbose_(verbose)
{}


Foam::topoSetSource::topoSetSource
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh)
{
    verbose(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::topoSetSource::verbose(const dictionary& dict)
{
    bool flag(verbose_);

    if (dict.readIfPresent("verbose", flag))
    {
        verbose_ = flag;
    }
}


// ************************************************************************* //

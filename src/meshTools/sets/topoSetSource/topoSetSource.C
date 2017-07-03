/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "topoSetSource.H"
#include "polyMesh.H"
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
Foam::topoSetSource::actionNames_
{
    { setAction::CLEAR, "clear" },
    { setAction::NEW, "new" },
    { setAction::INVERT, "invert" },
    { setAction::ADD, "add" },
    { setAction::DELETE, "delete" },
    { setAction::SUBSET, "subset" },
    { setAction::LIST, "list" },
    { setAction::REMOVE, "remove" },
};


const Foam::string Foam::topoSetSource::illegalSource_
(
    "Illegal topoSetSource name"
);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::topoSetSource> Foam::topoSetSource::New
(
    const word& topoSetSourceType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    auto cstrIter = wordConstructorTablePtr_->cfind(topoSetSourceType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown topoSetSource type "
            << topoSetSourceType << nl << nl
            << "Valid topoSetSource types :" << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<topoSetSource>(cstrIter()(mesh, dict));
}


Foam::autoPtr<Foam::topoSetSource> Foam::topoSetSource::New
(
    const word& topoSetSourceType,
    const polyMesh& mesh,
    Istream& is
)
{
    auto cstrIter = istreamConstructorTablePtr_->cfind(topoSetSourceType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown topoSetSource type "
            << topoSetSourceType << nl << nl
            << "Valid topoSetSource types :" << endl
            << istreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<topoSetSource>(cstrIter()(mesh, is));
}


Foam::Istream& Foam::topoSetSource::checkIs(Istream& is)
{
    if (is.good() && !is.eof())
    {
        return is;
    }
    else
    {
        FatalErrorInFunction
            << exit(FatalError);

        return is;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::topoSetSource::addOrDelete
(
    topoSet& set,
    const label celli,
    const bool add
) const
{
    if (add)
    {
        set.insert(celli);
    }
    else
    {
        set.erase(celli);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoSetSource::topoSetSource(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::topoSetSource::~topoSetSource()
{}


// ************************************************************************* //

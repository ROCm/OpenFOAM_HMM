/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "topoSetPointSource.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topoSetPointSource, 0);
    defineRunTimeSelectionTable(topoSetPointSource, word);
    defineRunTimeSelectionTable(topoSetPointSource, istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoSetPointSource::topoSetPointSource(const polyMesh& mesh)
:
    topoSetSource(mesh)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::topoSetPointSource> Foam::topoSetPointSource::New
(
    const word& sourceType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    auto cstrIter = wordConstructorTablePtr_->cfind(sourceType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown topoSetPointSource type "
            << sourceType << nl << nl
            << "Valid types :" << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<topoSetPointSource>(cstrIter()(mesh, dict));
}


Foam::autoPtr<Foam::topoSetPointSource> Foam::topoSetPointSource::New
(
    const word& sourceType,
    const polyMesh& mesh,
    Istream& is
)
{
    auto cstrIter = istreamConstructorTablePtr_->cfind(sourceType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown topoSetPointSource type "
            << sourceType << nl << nl
            << "Valid types :" << endl
            << istreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<topoSetPointSource>(cstrIter()(mesh, is));
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "token.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    typedef token::compound tokenCompound;
    defineTypeNameAndDebug(tokenCompound, 0);
    defineRunTimeSelectionTable(tokenCompound, Istream);
}

const char* const Foam::token::typeName = "token";
const Foam::token Foam::token::undefinedToken;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::token::parseError(const char* expected) const
{
    FatalIOError
        << "Parse error, expected a " << expected
        << ", found \n    " << info() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::token::compound::~compound()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::token::compound> Foam::token::compound::New
(
    const word& compoundType,
    Istream& is
)
{
    auto cstrIter = IstreamConstructorTablePtr_->cfind(compoundType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(is)
            << "Unknown compound type "
            << compoundType << nl << nl
            << "Valid compound types:" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << abort(FatalIOError);
    }

    return autoPtr<Foam::token::compound>(cstrIter()(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::token::compound::isCompound(const word& name)
{
    return
    (
        IstreamConstructorTablePtr_
     && IstreamConstructorTablePtr_->found(name)
    );
}


Foam::token::compound& Foam::token::transferCompoundToken(const Istream& is)
{
    if (type_ == tokenType::COMPOUND)
    {
        if (data_.compoundPtr->empty())
        {
            FatalIOErrorInFunction(is)
                << "compound has already been transferred from token\n    "
                << info() << abort(FatalIOError);
        }
        else
        {
            data_.compoundPtr->empty() = true;
        }

        return *data_.compoundPtr;
    }

    parseError("compound");
    return *data_.compoundPtr;
}


// ************************************************************************* //

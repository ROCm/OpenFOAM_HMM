/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

const Foam::token Foam::token::undefinedToken;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::token::parseError(const char* expected) const
{
    FatalIOError
        << "Parse error, expected a " << expected
        << ", found \n    " << info() << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::token::compound> Foam::token::compound::New
(
    const word& compoundType,
    Istream& is
)
{
    auto* ctorPtr = IstreamConstructorTable(compoundType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            is,
            "compound",
            compoundType,
            *IstreamConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    return autoPtr<Foam::token::compound>(ctorPtr(is));
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


Foam::token::compound& Foam::token::transferCompoundToken()
{
    if (type_ != tokenType::COMPOUND)
    {
        parseError("compound");
    }

    if (data_.compoundPtr->moved())
    {
        FatalErrorInFunction
            << "compound has already been transferred from token\n    "
            << info() << abort(FatalError);
    }
    else
    {
        data_.compoundPtr->moved(true);
    }

    return *data_.compoundPtr;
}


Foam::token::compound& Foam::token::transferCompoundToken(const Istream& is)
{
    if (type_ != tokenType::COMPOUND)
    {
        parseError("compound");
    }

    if (data_.compoundPtr->moved())
    {
        FatalIOErrorInFunction(is)
            << "compound has already been transferred from token\n    "
            << info() << abort(FatalIOError);
    }
    else
    {
        data_.compoundPtr->moved(true);
    }

    return *data_.compoundPtr;
}


// ************************************************************************* //

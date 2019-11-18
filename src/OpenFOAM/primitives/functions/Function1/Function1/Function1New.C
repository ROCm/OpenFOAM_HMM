/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "Constant.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::Function1<Type>> Foam::Function1<Type>::New
(
    const word& entryName,
    const dictionary& dict,
    const word& redirectType
)
{
    word modelType(redirectType);

    const entry* eptr = dict.findEntry(entryName, keyType::LITERAL);

    if (!eptr)
    {
        if (modelType.empty())
        {
            FatalIOErrorInFunction(dict)
                << "No Function1 dictionary entry: "
                << entryName << nl << nl
                << exit(FatalIOError);
        }
    }
    else if (eptr->isDict())
    {
        const dictionary& coeffsDict = eptr->dict();

        coeffsDict.readEntry
        (
            "type",
            modelType,
            keyType::LITERAL,
            redirectType.empty()  // mandatory when redirectType is empty
        );

        auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown Function1 type "
                << modelType << " for " << entryName
                << "\n\nValid Function1 types :\n"
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalIOError);
        }

        return cstrIter()(entryName, coeffsDict);
    }
    else
    {
        Istream& is = eptr->stream();

        token firstToken(is);

        if (!firstToken.isWord())
        {
            is.putBack(firstToken);
            return autoPtr<Function1<Type>>
            (
                new Function1Types::Constant<Type>(entryName, is)
            );
        }

        modelType = firstToken.wordToken();
    }


    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown Function1 type "
            << modelType << " for " << entryName
            << "\n\nValid Function1 types :\n"
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return cstrIter()
    (
        entryName,
        dict.optionalSubDict(entryName + "Coeffs")
    );
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    if (dict.isDict(entryName))
    {
        const dictionary& coeffsDict(dict.subDict(entryName));

        const word Function1Type
        (
            redirectType.empty()
          ? coeffsDict.get<word>("type")
          : coeffsDict.lookupOrDefault<word>("type", redirectType)
        );

        auto cstrIter = dictionaryConstructorTablePtr_->cfind(Function1Type);

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown Function1 type "
                << Function1Type << " for Function1 "
                << entryName << nl << nl
                << "Valid Function1 types :" << nl
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalIOError);
        }

        return cstrIter()(entryName, coeffsDict);
    }
    else
    {
        const dictionary::const_searcher finder
        (
            dict.csearch(entryName, keyType::REGEX)
        );

        word Function1Type;
        if (finder.found())
        {
            Istream& is = finder.ref().stream();

            token firstToken(is);

            if (!firstToken.isWord())
            {
                is.putBack(firstToken);
                return autoPtr<Function1<Type>>
                (
                    new Function1Types::Constant<Type>(entryName, is)
                );
            }
            else
            {
                Function1Type = firstToken.wordToken();
            }
        }
        else if (redirectType != word::null)
        {
            Function1Type = redirectType;
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find specification for Function1 "
                << entryName << nl << nl
                << "Valid Function1 types :" << nl
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalIOError);
        }

        auto cstrIter = dictionaryConstructorTablePtr_->cfind(Function1Type);

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown Function1 type "
                << Function1Type << " for Function1 "
                << entryName << nl << nl
                << "Valid Function1 types :" << nl
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalIOError);
        }

        return cstrIter()
        (
            entryName,
            dict.found(entryName + "Coeffs")
          ? dict.subDict(entryName + "Coeffs")
          : dict
        );
    }
}


// ************************************************************************* //

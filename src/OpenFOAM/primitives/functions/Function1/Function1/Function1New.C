/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "Constant.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::New
(
    const word& entryName,
    const entry* eptr,
    const dictionary& dict,
    const word& redirectType,
    const bool mandatory
)
{
    word modelType(redirectType);

    const dictionary* coeffs = (eptr ? eptr->dictPtr() : nullptr);

    if (coeffs)
    {
        // Dictionary entry

        DebugInFunction
            << "For " << entryName << " with dictionary entries: "
            << flatOutput(coeffs->toc()) << nl;

        coeffs->readEntry
        (
            "type",
            modelType,
            keyType::LITERAL,
            modelType.empty()  // "type" entry is mandatory if no 'redirect'
        );

        // Fallthrough
    }
    else if (eptr)
    {
        // Primitive entry
        // - word : the modelType
        // - non-word : value for constant function

        DebugInFunction
            << "For " << entryName << " with primitive entry" << nl;

        ITstream& is = eptr->stream();

        if (is.peek().isWord())
        {
            modelType = is.peek().wordToken();
        }
        else
        {
            // A value - compatibility for reading constant

            const Type constValue = pTraits<Type>(is);

            return autoPtr<Function1<Type>>
            (
                new Function1Types::Constant<Type>(entryName, constValue)
            );
        }

        // Fallthrough
    }


    if (modelType.empty())
    {
        // Entry missing

        if (mandatory)
        {
            FatalIOErrorInFunction(dict)
                << "Missing or invalid Function1 entry: "
                << entryName << nl
                << exit(FatalIOError);
        }

        return nullptr;
    }
    else if (!coeffs)
    {
        // Primitive entry. Coeffs dictionary is optional.

        const word& kw =
        (
            eptr
          ? eptr->keyword()  // Could be a compatibility lookup
          : entryName
        );

        coeffs = &dict.optionalSubDict(kw + "Coeffs", keyType::LITERAL);
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

    return cstrIter()(entryName, *coeffs);
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::New
(
    const word& entryName,
    const dictionary& dict,
    const word& redirectType,
    const bool mandatory
)
{
    return Function1<Type>::New
    (
        entryName,
        dict.findEntry(entryName, keyType::LITERAL),
        dict,
        redirectType,
        mandatory
    );
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::NewCompat
(
    const word& entryName,
    std::initializer_list<std::pair<const char*,int>> compat,
    const dictionary& dict,
    const word& redirectType,
    const bool mandatory
)
{
    return Function1<Type>::New
    (
        entryName,
        dict.findCompat(entryName, compat, keyType::LITERAL),
        dict,
        redirectType,
        mandatory
    );
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::New
(
    const word& entryName,
    const dictionary& dict,
    const bool mandatory
)
{
    return Function1<Type>::New(entryName, dict, word::null, mandatory);
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::NewIfPresent
(
    const word& entryName,
    const dictionary& dict,
    const word& redirectType
)
{
    // mandatory = false
    return Function1<Type>::New(entryName, dict, redirectType, false);
}


// ************************************************************************* //

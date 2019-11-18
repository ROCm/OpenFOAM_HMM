/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "ConstantField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::PatchFunction1<Type>> Foam::PatchFunction1<Type>::New
(
    const polyPatch& pp,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
{
    const entry* eptr = dict.findEntry(entryName, keyType::LITERAL);

    if (!eptr)
    {
        FatalIOErrorInFunction(dict)
            << "No PatchFunction1 dictionary entry: "
            << entryName << nl << nl
            << exit(FatalIOError);

        return nullptr;
    }
    else if (eptr->isDict())
    {
        const dictionary& coeffsDict = eptr->dict();

        const word modelType(coeffsDict.getWord("type", keyType::LITERAL));

        auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction(coeffsDict)
                << "Unknown PatchFunction1 type "
                << modelType << " for " << entryName
                << "\n\nValid PatchFunction1 types :\n"
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalIOError);
        }

        return cstrIter()(pp, modelType, entryName, coeffsDict, faceValues);
    }
    else
    {
        ITstream& is = eptr->stream();

        token firstToken(is);

        if (!firstToken.isWord())
        {
            // Backwards-compatibility for reading straight fields
            is.putBack(firstToken);

            const Type uniformValue = pTraits<Type>(is);
            const Field<Type> value
            (
                (faceValues ? pp.size() : pp.nPoints()),
                uniformValue
            );

            return autoPtr<PatchFunction1<Type>>
            (
                new PatchFunction1Types::ConstantField<Type>
                (
                    pp,
                    entryName,
                    true,           // uniform
                    uniformValue,   // uniform value
                    value,          // Supply value
                    dict,
                    faceValues
                )
            );
        }

        const word modelType = firstToken.wordToken();

        // Check to see if this looks like a normal field entry
        if (modelType == "uniform" || modelType == "nonuniform")
        {
            return autoPtr<PatchFunction1<Type>>
            (
                new PatchFunction1Types::ConstantField<Type>
                (
                    pp,
                    PatchFunction1Types::ConstantField<Type>::typeName,
                    entryName,
                    dict,
                    faceValues
                )
            );
        }


        auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown PatchFunction1 type "
                << modelType << " for " << entryName
                << "\n\nValid PatchFunction1 types :\n"
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalIOError);
        }

        return cstrIter()
        (
            pp,
            modelType,
            entryName,
            dict.optionalSubDict(entryName + "Coeffs"),
            faceValues
        );
    }
}


// ************************************************************************* //

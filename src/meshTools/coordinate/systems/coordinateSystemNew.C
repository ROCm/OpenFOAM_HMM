/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "objectRegistry.H"
#include "cartesianCS.H"
#include "indirectCS.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


//- Handle a 'coordinateSystem' sub-dictionary
// In 1806 and earlier, this was handled (rather poorly) in the
// coordinateSystem constructor itself.
const Foam::dictionary* Foam::coordinateSystem::subDictCompat
(
    const dictionary* dictPtr
)
{
    if (dictPtr)
    {
        // Non-recursive, no pattern matching in the search
        const auto finder =
            dictPtr->csearch(coordinateSystem::typeName_(), keyType::LITERAL);

        if (finder.isDict())
        {
            return finder.dictPtr();
        }
        else if (finder.found())
        {
            const word csName(finder.ref().stream());

            // Deprecated, unsupported syntax
            if (error::master())
            {
                std::cerr
                    << "--> FOAM IOWarning :" << nl
                    << "    Ignoring 'coordinateSystem' as a keyword."
                    " Perhaps you meant this instead?" << nl
                    << '{' << nl
                    << "    type " << coordSystem::indirect::typeName_()
                    << ';' << nl
                    << "    name " << csName << ';' << nl
                    << '}' << nl
                    << std::endl;

                error::warnAboutAge("syntax change", 1806);
            }
        }
    }

    return dictPtr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    word modelType,
    const objectRegistry& obr,
    const dictionary& dict
)
{
    if (modelType.empty())
    {
        modelType = coordSystem::cartesian::typeName_();
    }

    {
        auto* ctorPtr = registryConstructorTable(modelType);
        if (ctorPtr)
        {
            return autoPtr<coordinateSystem>(ctorPtr(obr, dict));
        }
    }

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    // Everything with a registry constructor also has a dictionary
    // constructor, so just need to print those.
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "coordinate system",
             modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<coordinateSystem>(ctorPtr(dict));
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    word modelType,
    const dictionary& dict
)
{
    if (modelType.empty())
    {
        modelType = coordSystem::cartesian::typeName_();
    }

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "coordinate system",
             modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<coordinateSystem>(ctorPtr(dict));
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& dictName
)
{
    const dictionary* dictPtr = &dict;

    if (dictName.size())
    {
        dictPtr = &(dictPtr->subDict(dictName));
    }
    else
    {
        // Use 'coordinateSystem' subDict if present
        dictPtr = coordinateSystem::subDictCompat(dictPtr);
    }

    word modelType = dictPtr->getOrDefault<word>
    (
        "type",
        coordSystem::cartesian::typeName_()
    );

    return coordinateSystem::New(modelType, obr, *dictPtr);
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    const dictionary& dict,
    const word& dictName
)
{
    const dictionary* dictPtr = &dict;

    if (dictName.size())
    {
        dictPtr = &(dictPtr->subDict(dictName));
    }
    else
    {
        // Use 'coordinateSystem' subDict if present
        dictPtr = coordinateSystem::subDictCompat(dictPtr);
    }

    const word modelType = dictPtr->getOrDefault<word>
    (
        "type",
        coordSystem::cartesian::typeName_()
    );

    return coordinateSystem::New(modelType, *dictPtr);
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New(Istream& is)
{
    const word csName(is);
    const dictionary dict(is);

    auto cs = coordinateSystem::New(dict, word::null);
    cs->rename(csName);

    return cs;
}


// ************************************************************************* //

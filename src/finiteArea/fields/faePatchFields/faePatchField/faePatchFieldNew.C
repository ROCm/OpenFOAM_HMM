/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::faePatchField<Type>> Foam::faePatchField<Type>::New
(
    const word& patchFieldType,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF
)
{
    DebugInFunction << "Constructing faePatchField" << endl;

    auto* ctorPtr = patchConstructorTable(patchFieldType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "patchField",
            patchFieldType,
            *patchConstructorTablePtr_
        ) << exit(FatalError);
    }

    auto* patchTypeCtor = patchConstructorTable(p.type());

    if (patchTypeCtor)
    {
        return patchTypeCtor(p, iF);
    }
    else
    {
        return ctorPtr(p, iF);
    }
}


template<class Type>
Foam::tmp<Foam::faePatchField<Type>> Foam::faePatchField<Type>::New
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const dictionary& dict
)
{
    DebugInFunction << "Constructing faePatchField" << endl;

    const word patchFieldType(dict.get<word>("type"));

    auto* ctorPtr = dictionaryConstructorTable(patchFieldType);

    if (!ctorPtr)
    {
        if (!disallowGenericFaePatchField)
        {
            ctorPtr = dictionaryConstructorTable("generic");
        }

        if (!ctorPtr)
        {
            FatalIOErrorInFunction(dict)
                << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << nl << nl
                << "Valid patchField types are :" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    auto* patchTypeCtor = dictionaryConstructorTable(p.type());

    if (patchTypeCtor && patchTypeCtor != ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "inconsistent patch and patchField types for \n"
            << "    patch type " << p.type()
            << " and patchField type " << patchFieldType
            << exit(FatalIOError);
    }

    return ctorPtr(p, iF, dict);
}


template<class Type>
Foam::tmp<Foam::faePatchField<Type>> Foam::faePatchField<Type>::New
(
    const faePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper& pfMapper
)
{
    DebugInFunction << "Constructing faePatchField<Type>" << endl;

    auto* ctorPtr = patchMapperConstructorTable(ptf.type());

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "patchField",
            ptf.type(),
            *patchMapperConstructorTablePtr_
        ) << exit(FatalError);
    }

    auto* patchTypeCtor = patchMapperConstructorTable(p.type());

    if (patchTypeCtor)
    {
        return patchTypeCtor(ptf, p, iF, pfMapper);
    }

    return ctorPtr(ptf, p, iF, pfMapper);
}


// ************************************************************************* //

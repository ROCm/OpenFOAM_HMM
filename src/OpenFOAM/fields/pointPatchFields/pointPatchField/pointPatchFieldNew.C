/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
{
    DebugInFunction << "Constructing pointPatchField<Type>" << endl;

    auto* ctorPtr = pointPatchConstructorTable(patchFieldType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "patchFieldType",
            patchFieldType,
            *pointPatchConstructorTablePtr_
        ) << exit(FatalError);
    }

    autoPtr<pointPatchField<Type>> pfPtr(ctorPtr(p, iF));

    if
    (
        actualPatchType.empty()
     || actualPatchType != p.type()
    )
    {
        if (pfPtr().constraintType() != p.constraintType())
        {
            // Incompatible (constraint-wise) with the patch type
            // - use default constraint type

            auto* patchTypeCtor = pointPatchConstructorTable(p.type());

            if (!patchTypeCtor)
            {
                FatalErrorInFunction
                    << "Inconsistent patch and patchField types for\n"
                    << "    patch type " << p.type()
                    << " and patchField type " << patchFieldType
                    << exit(FatalError);
            }

            return patchTypeCtor(p, iF);
        }
    }
    else
    {
        if (pointPatchConstructorTablePtr_->found(p.type()))
        {
            pfPtr().patchType() = actualPatchType;
        }
    }

    return pfPtr;
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const word& patchFieldType,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
{
    DebugInFunction << "Constructing pointPatchField<Type>" << endl;

    const word patchFieldType(dict.get<word>("type"));

    auto* ctorPtr = dictionaryConstructorTable(patchFieldType);

    if (!ctorPtr)
    {
        if (!disallowGenericPointPatchField)
        {
            ctorPtr = dictionaryConstructorTable("generic");
        }

        if (!ctorPtr)
        {
            FatalIOErrorInFunction(dict)
                << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << nl << nl
                << "Valid patchField types :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    // Construct (but not necessarily returned)
    autoPtr<pointPatchField<Type>> pfPtr(ctorPtr(p, iF, dict));

    if
    (
        !dict.found("patchType")
     || dict.get<word>("patchType") != p.type()
    )
    {
        if (pfPtr().constraintType() != p.constraintType())
        {
            // Incompatible (constraint-wise) with the patch type
            // - use default constraint type

            auto* patchTypeCtor = dictionaryConstructorTable(p.type());

            if (!patchTypeCtor)
            {
                FatalIOErrorInFunction(dict)
                    << "Inconsistent patch and patchField types for\n"
                    << "    patch type " << p.type()
                    << " and patchField type " << patchFieldType
                    << exit(FatalIOError);
            }

            return patchTypeCtor(p, iF, dict);
        }
    }

    return pfPtr;
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const pointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& pfMapper
)
{
    DebugInFunction << "Constructing pointPatchField<Type>" << endl;

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

    return ctorPtr(ptf, p, iF, pfMapper);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
Foam::tmp<Foam::faPatchField<Type>> Foam::faPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
{
    DebugInFunction
        << "Constructing faPatchField<Type> "
        << "patchFieldType:" << patchFieldType
        << "actualPatchType:" << actualPatchType
        << "p.Type():" << p.type()
        << endl;

    auto cstrIter = patchConstructorTablePtr_->cfind(patchFieldType);

    if (!cstrIter.found())
    {
        FatalErrorInLookup
        (
            "patchField",
            patchFieldType,
            *patchConstructorTablePtr_
        ) << exit(FatalError);
    }

    auto patchTypeCstrIter = patchConstructorTablePtr_->cfind(p.type());

    if
    (
        actualPatchType == word::null
     || actualPatchType != p.type()
    )
    {
        if (patchTypeCstrIter.found())
        {
            return patchTypeCstrIter()(p, iF);
        }
        else
        {
            return cstrIter()(p, iF);
        }
    }


    tmp<faPatchField<Type>> tfap = cstrIter()(p, iF);

    // Check if constraint type override and store patchType if so
    if (patchTypeCstrIter.found())
    {
        tfap.ref().patchType() = actualPatchType;
    }
    return tfap;
}


template<class Type>
Foam::tmp<Foam::faPatchField<Type>> Foam::faPatchField<Type>::New
(
    const word& patchFieldType,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}


template<class Type>
Foam::tmp<Foam::faPatchField<Type>> Foam::faPatchField<Type>::New
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
{
    DebugInFunction << "Constructing faPatchField<Type>" << endl;

    const word patchFieldType(dict.get<word>("type"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(patchFieldType);

    if (!cstrIter.found())
    {
        if (!disallowGenericFaPatchField)
        {
            cstrIter = dictionaryConstructorTablePtr_->cfind("generic");
        }

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << nl << nl
                << "Valid patchField types are :" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    auto patchTypeCstrIter = dictionaryConstructorTablePtr_->cfind(p.type());

    if (patchTypeCstrIter.found() && *patchTypeCstrIter != *cstrIter)
    {
        FatalIOErrorInFunction(dict)
            << "inconsistent patch and patchField types for \n"
            << "    patch type " << p.type()
            << " and patchField type " << patchFieldType
            << exit(FatalIOError);
    }

    return cstrIter()(p, iF, dict);
}


template<class Type>
Foam::tmp<Foam::faPatchField<Type>> Foam::faPatchField<Type>::New
(
    const faPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& pfMapper
)
{
    DebugInFunction << "Constructing faPatchField<Type>" << endl;

    auto cstrIter = patchMapperConstructorTablePtr_->cfind(ptf.type());

    if (!cstrIter.found())
    {
        FatalErrorInLookup
        (
            "patchField",
            ptf.type(),
            *patchMapperConstructorTablePtr_
        ) << exit(FatalError);
    }

    auto patchTypeCstrIter = patchMapperConstructorTablePtr_->cfind(p.type());

    if (patchTypeCstrIter.found())
    {
        return patchTypeCstrIter()(ptf, p, iF, pfMapper);
    }

    return cstrIter()(ptf, p, iF, pfMapper);
}


// ************************************************************************* //

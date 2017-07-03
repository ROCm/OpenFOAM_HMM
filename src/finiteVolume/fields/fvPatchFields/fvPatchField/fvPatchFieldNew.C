/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    if (debug)
    {
        InfoInFunction
            << "patchFieldType = " << patchFieldType
            << " : " << p.type()
            << endl;
    }

    auto cstrIter = patchConstructorTablePtr_->cfind(patchFieldType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown patchField type "
            << patchFieldType << nl << nl
            << "Valid patchField types :" << endl
            << patchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
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
    else
    {
        tmp<fvPatchField<Type>> tfvp = cstrIter()(p, iF);

        // Check if constraint type override and store patchType if so
        if (patchTypeCstrIter.found())
        {
            tfvp.ref().patchType() = actualPatchType;
        }
        return tfvp;
    }
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const word& patchFieldType,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
{
    const word patchFieldType(dict.lookup("type"));

    if (debug)
    {
        InfoInFunction
            << "patchFieldType = " << patchFieldType
            << endl;
    }

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(patchFieldType);

    if (!cstrIter.found())
    {
        if (!disallowGenericFvPatchField)
        {
            cstrIter = dictionaryConstructorTablePtr_->cfind("generic");
        }

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << nl << nl
                << "Valid patchField types :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    if
    (
       !dict.found("patchType")
     || word(dict.lookup("patchType")) != p.type()
    )
    {
        auto patchTypeCstrIter
            = dictionaryConstructorTablePtr_->cfind(p.type());

        if (patchTypeCstrIter.found() && patchTypeCstrIter() != cstrIter())
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "inconsistent patch and patchField types for \n"
                   "    patch type " << p.type()
                << " and patchField type " << patchFieldType
                << exit(FatalIOError);
        }
    }

    return cstrIter()(p, iF, dict);
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::New
(
    const fvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvPatchField<Type>" << endl;
    }

    auto cstrIter = patchMapperConstructorTablePtr_->cfind(ptf.type());

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown patchField type " << ptf.type() << nl << nl
            << "Valid patchField types :" << endl
            << patchMapperConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(ptf, p, iF, pfMapper);
}


// ************************************************************************* //

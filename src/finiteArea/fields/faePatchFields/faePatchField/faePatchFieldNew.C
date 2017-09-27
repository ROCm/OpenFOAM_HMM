/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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
    DebugInFunction
        << "constructing faePatchField<Type>"
        << endl;

    auto cstrIter = patchConstructorTablePtr_->cfind(patchFieldType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown patchTypefield type " << patchFieldType << nl << nl
            << "Valid patchField types are :" << nl
            << patchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    auto patchTypeCstrIter = patchConstructorTablePtr_->cfind(p.type());

    if (patchTypeCstrIter.found())
    {
        return patchTypeCstrIter()(p, iF);
    }
    else
    {
        return cstrIter()(p, iF);
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
    DebugInFunction
        << "constructing faePatchField<Type>"
        << endl;

    word patchFieldType(dict.lookup("type"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(patchFieldType);

    if (!cstrIter.found())
    {
        if (!disallowGenericFaePatchField)
        {
            cstrIter = dictionaryConstructorTablePtr_->find("generic");
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
Foam::tmp<Foam::faePatchField<Type>> Foam::faePatchField<Type>::New
(
    const faePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper& pfMapper
)
{
    DebugInFunction
        << "constructing faePatchField<Type>"
        << endl;

    auto cstrIter = patchMapperConstructorTablePtr_->cfind(ptf.type());

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "unknown patchTypefield type " << ptf.type() << endl << endl
            << "Valid patchField types are :" << endl
            << patchMapperConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    typename patchMapperConstructorTable::iterator
        patchTypeCstrIter = patchMapperConstructorTablePtr_->find(p.type());

    if (patchTypeCstrIter.found())
    {
        return patchTypeCstrIter()(ptf, p, iF, pfMapper);
    }
    else
    {
        return cstrIter()(ptf, p, iF, pfMapper);
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "calculatedPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::pointPatchField<Type>::calculatedType()
{
    return calculatedPointPatchField<Type>::typeName;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::calculatedPointPatchField<Type>::calculatedPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    pointPatchField<Type>(p, iF)
{}


template<class Type>
Foam::calculatedPointPatchField<Type>::calculatedPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    pointPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::calculatedPointPatchField<Type>::calculatedPointPatchField
(
    const calculatedPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    pointPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::calculatedPointPatchField<Type>::calculatedPointPatchField
(
    const calculatedPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    pointPatchField<Type>(ptf, iF)
{}


template<class Type>
template<class Type2>
Foam::autoPtr<Foam::pointPatchField<Type>>
Foam::pointPatchField<Type>::NewCalculatedType
(
    const pointPatchField<Type2>& pf
)
{
    auto* patchTypeCtor = pointPatchConstructorTable(pf.patch().type());

    if (patchTypeCtor)
    {
        return autoPtr<pointPatchField<Type>>
        (
            patchTypeCtor
            (
                pf.patch(),
                Field<Type>::null()
            )
        );
    }
    else
    {
        return autoPtr<pointPatchField<Type>>
        (
            new calculatedPointPatchField<Type>
            (
                pf.patch(),
                Field<Type>::null()
            )
        );
    }
}


// ************************************************************************* //

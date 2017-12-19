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

#include "calculatedFaPatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::faPatchField<Type>::calculatedType()
{
    return calculatedFaPatchField<Type>::typeName;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::calculatedFaPatchField<Type>::calculatedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF)
{}


template<class Type>
Foam::calculatedFaPatchField<Type>::calculatedFaPatchField
(
    const calculatedFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::calculatedFaPatchField<Type>::calculatedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, Field<Type>("value", dict, p.size()))
{}


template<class Type>
Foam::calculatedFaPatchField<Type>::calculatedFaPatchField
(
    const calculatedFaPatchField<Type>& ptf
)
:
    faPatchField<Type>(ptf)
{}


template<class Type>
Foam::calculatedFaPatchField<Type>::calculatedFaPatchField
(
    const calculatedFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf, iF)
{}


template<class Type>
template<class Type2>
Foam::tmp<Foam::faPatchField<Type>> Foam::faPatchField<Type>::NewCalculatedType
(
    const faPatchField<Type2>& pf
)
{
    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTablePtr_->find(pf.patch().type());

    if (patchTypeCstrIter.found())
    {
        return patchTypeCstrIter()
        (
            pf.patch(),
            DimensionedField<Type, areaMesh>::null()
        );
    }
    else
    {
        return tmp<faPatchField<Type>>
        (
            new calculatedFaPatchField<Type>
            (
                pf.patch(),
                DimensionedField<Type, areaMesh>::null()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorInFunction
        << "valueInternalCoeffs cannot be called for a calculatedFaPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << exit(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorInFunction
        << "valueBoundaryCoeffs cannot be called for a calculatedFaPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << exit(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFaPatchField<Type>::gradientInternalCoeffs() const
{
    FatalErrorInFunction
        << "gradientInternalCoeffs cannot be called for a "
           "calculatedFaPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << exit(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    FatalErrorInFunction
        << "\n    "
           "gradientBoundaryCoeffs cannot be called for a "
           "calculatedFaPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << exit(FatalError);

    return *this;
}


template<class Type>
void Foam::calculatedFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //

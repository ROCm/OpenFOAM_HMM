/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "uniformFixedValueFaPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::uniformFixedValueFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(p, iF),
    uniformValue_(nullptr)
{}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::uniformFixedValueFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFaPatchField<Type>(p, iF),
    uniformValue_(Function1<Type>::New("uniformValue", dict))
{
   if (dict.found("value"))
   {
       faPatchField<Type>::operator==
       (
           Field<Type>("value", dict, p.size())
       );
   }
   else
   {
       this->evaluate();
   }
}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::uniformFixedValueFaPatchField
(
    const uniformFixedValueFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedValueFaPatchField<Type>(ptf, p, iF, mapper),
    uniformValue_(ptf.uniformValue_.clone())
{}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::uniformFixedValueFaPatchField
(
    const uniformFixedValueFaPatchField<Type>& ptf
)
:
    fixedValueFaPatchField<Type>(ptf),
    uniformValue_(ptf.uniformValue_.clone())
{}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::uniformFixedValueFaPatchField
(
    const uniformFixedValueFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedValueFaPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    faPatchField<Type>::operator==(uniformValue_->value(t));
    fixedValueFaPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformFixedValueFaPatchField<Type>::write
(
    Ostream& os
) const
{
    faPatchField<Type>::write(os);
    uniformValue_->writeData(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //

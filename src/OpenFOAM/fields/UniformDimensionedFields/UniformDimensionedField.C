/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "UniformDimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const IOobject& io,
    const dimensioned<Type>& dt
)
:
    regIOobject(io),
    dimensioned<Type>(dt)
{
    // Read value
    readHeaderOk(IOstream::BINARY, typeName);
}


template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const UniformDimensionedField<Type>& rdt
)
:
    regIOobject(rdt),
    dimensioned<Type>(rdt)
{}


template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const IOobject& io
)
:
    regIOobject(io),
    dimensioned<Type>(regIOobject::name(), dimless, Zero)
{
    // For if MUST_READ_IF_MODIFIED
    addWatch();

    // Read unless NO_READ
    readHeaderOk(IOstream::BINARY, typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::UniformDimensionedField<Type>::~UniformDimensionedField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::UniformDimensionedField<Type>::readData(Istream& is)
{
    dictionary dict(is);
    scalar multiplier(1);
    this->dimensions().read
    (
        dict.lookup("dimensions", keyType::LITERAL),
        multiplier
    );

    dict.readEntry("value", this->value());
    this->value() *= multiplier;

    return is.good();
}


template<class Type>
bool Foam::UniformDimensionedField<Type>::writeData(Ostream& os) const
{
    scalar multiplier(1);
    os.writeKeyword("dimensions");
    this->dimensions().write(os, multiplier) << token::END_STATEMENT << nl;
    os.writeEntry("value", this->value()/multiplier) << nl;

    return os.good();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::UniformDimensionedField<Type>::operator=
(
    const UniformDimensionedField<Type>& rhs
)
{
    dimensioned<Type>::operator=(rhs);
}


template<class Type>
void Foam::UniformDimensionedField<Type>::operator=
(
    const dimensioned<Type>& rhs
)
{
    dimensioned<Type>::operator=(rhs);
}


// ************************************************************************* //

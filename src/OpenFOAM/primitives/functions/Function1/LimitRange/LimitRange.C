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

#include "LimitRange.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::LimitRange<Type>::read(const dictionary& coeffs)
{
    min_ = coeffs.get<scalar>("min");
    max_ = coeffs.get<scalar>("max");
    value_ = Function1<Type>::New("value", coeffs);
}


template<class Type>
Foam::Function1Types::LimitRange<Type>::LimitRange
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName)
{
    read(dict);
}


template<class Type>
Foam::Function1Types::LimitRange<Type>::LimitRange(const LimitRange<Type>& rhs)
:
    Function1<Type>(rhs),
    min_(rhs.min_),
    max_(rhs.max_),
    value_(rhs.value_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::LimitRange<Type>::writeEntries(Ostream& os) const
{
    os.writeEntry("min", min_);
    os.writeEntry("max", max_);
    value_->writeData(os);
}


template<class Type>
void Foam::Function1Types::LimitRange<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os  << token::END_STATEMENT << nl;

    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //

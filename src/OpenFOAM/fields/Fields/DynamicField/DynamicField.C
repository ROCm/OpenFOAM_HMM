/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "DynamicField.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, int SizeMin>
Foam::DynamicField<T, SizeMin>::DynamicField(Istream& is)
:
    Field<T>(is),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
Foam::tmp<Foam::DynamicField<T, SizeMin>>
Foam::DynamicField<T, SizeMin>::clone() const
{
    return tmp<DynamicField<T, SizeMin>>
    (
        new DynamicField<T, SizeMin>(*this)
    );
}


// * * * * * * * * * * * * * * * IOstream Operator * * * * * * * * * * * * * //

template<class T, int SizeMin>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DynamicField<T, SizeMin>& lst
)
{
    os << static_cast<const Field<T>&>(lst);
    return os;
}


template<class T, int SizeMin>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    DynamicField<T, SizeMin>& lst
)
{
    is >> static_cast<Field<T>&>(lst);
    lst.capacity_ = lst.Field<T>::size();

    return is;
}


// ************************************************************************* //

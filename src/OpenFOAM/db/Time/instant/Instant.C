/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "Instant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::Instant<T>::Instant()
:
    val_(0),
    key_()
{}


template<class T>
Foam::Instant<T>::Instant::Instant(scalar val, const T& key)
:
    val_(val),
    key_(key)
{}


template<class T>
Foam::Instant<T>::Instant::Instant(scalar val, T&& key)
:
    val_(val),
    key_(std::move(key))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
bool Foam::Instant<T>::equal(scalar val) const
{
    return ((val_ > val - SMALL) && (val_ < val + SMALL));
}


template<class T>
template<class T2>
bool Foam::Instant<T>::equal(const Instant<T2>& other) const
{
    return this->equal(other.value());
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class T1, class T2>
bool Foam::operator==(const Instant<T1>& a, const Instant<T2>& b)
{
    return a.equal(b.value());
}


template<class T1, class T2>
bool Foam::operator!=(const Instant<T1>& a, const Instant<T2>& b)
{
    return !a.equal(b.value());
}


template<class T1, class T2>
bool Foam::operator<(const Instant<T1>& a, const Instant<T2>& b)
{
    return a.value() < b.value();
}


template<class T1, class T2>
bool Foam::operator>(const Instant<T1>& a, const Instant<T2>& b)
{
    return a.value() > b.value();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Foam::Istream& Foam::operator>>(Istream& is, Instant<T>& inst)
{
    is >> inst.value() >> inst.name();
    return is;
}


template<class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const Instant<T>& inst)
{
    os << inst.value() << '\t' << inst.name();
    return os;
}


// ************************************************************************* //

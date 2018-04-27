/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "clockValue.H"
#include "IOstreams.H"
#include <sys/time.h>

// * * * * * * * * * * * * * * * * Local Data  * * * * * * * * * * * * * * * //

namespace
{

    constexpr int factorMicro = (1000000);
    constexpr int factorMicro2 = (500000);

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clockValue::clockValue()
:
    clockValue(false)
{}


Foam::clockValue::clockValue(bool useNow)
{
    if (useNow)
    {
        update();
    }
    else
    {
        clear();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::clockValue::clear()
{
    value_.tv_sec  = 0;
    value_.tv_usec = 0;
}


void Foam::clockValue::update()
{
    gettimeofday(&value_, 0);
}


Foam::clockValue Foam::clockValue::elapsed() const
{
    return (now() -= *this);
}


long Foam::clockValue::seconds() const
{
    long sec = value_.tv_sec;
    if (sec > 0 && value_.tv_usec > factorMicro2)
    {
        ++sec;
    }

    return sec;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::clockValue::operator double () const
{
    return (value_.tv_sec + 1e-6*value_.tv_usec);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::clockValue& Foam::clockValue::operator-=(const clockValue& rhs)
{
    const value_type& b = rhs.value_;

    value_.tv_sec -= b.tv_sec;

    if (value_.tv_usec < b.tv_usec)
    {
        --(value_.tv_sec);
        value_.tv_usec += factorMicro;
    }

    value_.tv_usec -= b.tv_usec;

    return *this;
}


Foam::clockValue& Foam::clockValue::operator+=(const clockValue& rhs)
{
    const value_type& b = rhs.value_;

    // Microseconds first
    value_.tv_usec += b.tv_usec;

    value_.tv_sec  += b.tv_sec + (value_.tv_usec / factorMicro);
    value_.tv_usec %= factorMicro;

    return *this;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::clockValue Foam::operator-(const clockValue& a, const clockValue& b)
{
    clockValue result(a);
    result -= b;

    return result;
}


Foam::clockValue Foam::operator+(const clockValue& a, const clockValue& b)
{
    clockValue result(a);
    result += b;

    return result;
}


// ************************************************************************* //

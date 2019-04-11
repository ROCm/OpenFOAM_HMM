/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include <sstream>
#include <iomanip>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clockValue::clockValue()
:
    value_(value_type::zero())
{}


Foam::clockValue::clockValue(const value_type& value)
:
    value_(value)
{}


Foam::clockValue::clockValue(bool useNow)
:
    value_(value_type::zero())
{
    if (useNow)
    {
        update();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::clockValue::clear()
{
    value_ = value_type::zero();
}


void Foam::clockValue::update()
{
    value_ = std::chrono::high_resolution_clock::now().time_since_epoch();
}


Foam::clockValue Foam::clockValue::elapsed() const
{
    return (now() -= *this);
}


long Foam::clockValue::seconds() const
{
    return std::chrono::duration_cast<std::chrono::seconds>(value_).count();
}


std::string Foam::clockValue::str() const
{
    std::ostringstream os;

    // seconds
    const unsigned long ss =
         std::chrono::duration_cast<std::chrono::seconds>(value_).count();

    // days
    const auto dd = (ss / 86400);

    // hours
    const int hh = ((ss / 3600) % 24);

    if (dd) os << dd << '-';

    if (dd || hh)
    {
        os  << std::setw(2) << std::setfill('0')
            << hh << ':';
    }

    // minutes
    os  << std::setw(2) << std::setfill('0')
        << ((ss / 60) % 60) << ':';

    // seconds
    os  << std::setw(2) << std::setfill('0')
        << (ss % 60);

    // milliseconds. As none or 3 decimal places
    const long ms =
    (
        std::chrono::duration_cast<std::chrono::milliseconds>(value_).count()
      - (ss * 1000)
    );

    if (ms > 0)
    {
        os  << '.' << std::setw(3) << std::setfill('0') << ms;
    }

    return os.str();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::clockValue::operator double () const
{
    return
    (
        (double(value_.count()) * value_type::period::num)
      / value_type::period::den
    );
}


Foam::clockValue& Foam::clockValue::operator-=(const clockValue& rhs)
{
    value_ -= rhs.value_;
    return *this;
}


Foam::clockValue& Foam::clockValue::operator+=(const clockValue& rhs)
{
    value_ += rhs.value_;
    return *this;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::clockValue Foam::operator-(const clockValue& a, const clockValue& b)
{
    return clockValue(a.value() - b.value());
}


Foam::clockValue Foam::operator+(const clockValue& a, const clockValue& b)
{
    return clockValue(a.value() + b.value());
}


// ************************************************************************* //

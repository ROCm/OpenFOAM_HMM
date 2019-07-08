/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "int64.H"
#include "error.H"
#include "parsing.H"
#include "IOstreams.H"
#include <cinttypes>

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

int64_t Foam::readInt64(const char* buf)
{
    char *endptr = nullptr;
    errno = 0;
    const intmax_t parsed = ::strtoimax(buf, &endptr, 10);

    const int64_t val = int64_t(parsed);

    const parsing::errorType err =
    (
        (parsed < INT64_MIN || parsed > INT64_MAX)
      ? parsing::errorType::RANGE
      : parsing::checkConversion(buf, endptr)
    );

    if (err != parsing::errorType::NONE)
    {
        FatalIOErrorInFunction("unknown")
            << parsing::errorNames[err] << " '" << buf << "'"
            << exit(FatalIOError);
    }

    return val;
}


bool Foam::readInt64(const char* buf, int64_t& val)
{
    char *endptr = nullptr;
    errno = 0;
    const intmax_t parsed = ::strtoimax(buf, &endptr, 10);

    val = int64_t(parsed);

    return
    (
        (parsed < INT64_MIN || parsed > INT64_MAX)
      ? false
      : (parsing::checkConversion(buf, endptr) == parsing::errorType::NONE)
    );
}


Foam::Istream& Foam::operator>>(Istream& is, int64_t& val)
{
    token t(is);

    if (!t.good())
    {
        FatalIOErrorInFunction(is)
            << "Bad token - could not get int64"
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        val = int64_t(t.labelToken());
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong token type - expected label (int64), found "
            << t.info()
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


int64_t Foam::readInt64(Istream& is)
{
    int64_t val(0);
    is >> val;

    return val;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const int64_t val)
{
    os.write(label(val));
    os.check(FUNCTION_NAME);
    return os;
}


#if defined(__APPLE__)
Foam::Istream& Foam::operator>>(Istream& is, long& val)
{
    return operator>>(is, reinterpret_cast<int64_t&>(val));
}

Foam::Ostream& Foam::operator<<(Ostream& os, const long val)
{
    return (os << int64_t(val));
}
#endif


// ************************************************************************* //

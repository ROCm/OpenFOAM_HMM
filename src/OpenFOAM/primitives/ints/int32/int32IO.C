/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "int32.H"
#include "error.H"
#include "parsing.H"
#include "IOstreams.H"
#include <cinttypes>

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

int32_t Foam::readInt32(const char* buf)
{
    char *endptr = nullptr;
    errno = 0;
    const intmax_t parsed = ::strtoimax(buf, &endptr, 10);

    const int32_t val = int32_t(parsed);

    const parsing::errorType err =
    (
        (parsed < INT32_MIN || parsed > INT32_MAX)
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


bool Foam::readInt32(const char* buf, int32_t& val)
{
    char *endptr = nullptr;
    errno = 0;
    const intmax_t parsed = ::strtoimax(buf, &endptr, 10);

    val = int32_t(parsed);

    return
    (
        (parsed < INT32_MIN || parsed > INT32_MAX)
      ? false
      : (parsing::checkConversion(buf, endptr) == parsing::errorType::NONE)
    );
}


Foam::Istream& Foam::operator>>(Istream& is, int32_t& val)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        val = int32_t(t.labelToken());
    }
    else
    {
        is.setBad();
        FatalIOErrorInFunction(is)
            << "wrong token type - expected int32_t, found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


int32_t Foam::readInt32(Istream& is)
{
    int32_t val;
    is >> val;

    return val;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const int32_t val)
{
    os.write(label(val));
    os.check(FUNCTION_NAME);
    return os;
}


#if WM_ARCH_OPTION == 32
Foam::Istream& Foam::operator>>(Istream& is, long& val)
{
    return operator>>(is, reinterpret_cast<int32_t&>(val));
}

Foam::Ostream& Foam::operator<<(Ostream& os, const long val)
{
    os << int32_t(val);
    return os;
}
#endif


// ************************************************************************* //

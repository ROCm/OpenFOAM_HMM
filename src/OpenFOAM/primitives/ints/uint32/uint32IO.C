/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "uint32.H"
#include "parsing.H"
#include "IOstreams.H"
#include <cinttypes>

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

uint32_t Foam::readUint32(const char* buf)
{
    char *endptr = nullptr;
    errno = 0;
    const uintmax_t parsed = ::strtoumax(buf, &endptr, 10);

    const uint32_t val = uint32_t(parsed);

    const parsing::errorType err =
    (
        (parsed > UINT32_MAX)
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


bool Foam::readUint32(const char* buf, uint32_t& val)
{
    char *endptr = nullptr;
    errno = 0;
    const uintmax_t parsed = ::strtoumax(buf, &endptr, 10);

    val = uint32_t(parsed);

    return
    (
        (parsed > UINT32_MAX)
      ? false
      : (parsing::checkConversion(buf, endptr) == parsing::errorType::NONE)
    );
}


Foam::Istream& Foam::operator>>(Istream& is, uint32_t& val)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        val = uint32_t(t.labelToken());
    }
    else
    {
        is.setBad();
        FatalIOErrorInFunction(is)
            << "wrong token type - expected uint32_t, found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


uint32_t Foam::readUint32(Istream& is)
{
    uint32_t val;
    is >> val;

    return val;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const uint32_t val)
{
    os.write(label(val));
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //

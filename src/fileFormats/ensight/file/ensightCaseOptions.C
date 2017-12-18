/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "ensightCase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightCase::options::options(IOstream::streamFormat format)
:
    format_(format),
    width_(0),
    mask_(),
    printf_(),
    overwrite_(false),
    nodeValues_(false),
    separateCloud_(false)
{
    width(8); // ensures that the mask and printf-format are also resized
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::IOstream::streamFormat Foam::ensightCase::options::format() const
{
    return format_;
}

const Foam::word& Foam::ensightCase::options::mask() const
{
    return mask_;
}


Foam::word Foam::ensightCase::options::padded(const label i) const
{
    // As per Foam::name, but with fixed length
    char buf[32];

    ::snprintf(buf, 32, printf_.c_str(), static_cast<int>(i));
    buf[31] = 0;

    // no stripping required
    return word(buf, false);
}


Foam::label Foam::ensightCase::options::width() const
{
    return width_;
}


void Foam::ensightCase::options::width(const label n)
{
    // enforce min/max sanity limits
    if (n < 1 || n > 31)
    {
        return;
    }

    // set mask accordingly
    mask_.resize(n, '*');

    // appropriate printf format
    printf_ = "%0" + Foam::name(n) + "d";
}



bool Foam::ensightCase::options::overwrite() const
{
    return overwrite_;
}


void Foam::ensightCase::options::overwrite(bool b)
{
    overwrite_ = b;
}


bool Foam::ensightCase::options::nodeValues() const
{
    return nodeValues_;
}


void Foam::ensightCase::options::nodeValues(bool b)
{
    nodeValues_ = b;
}


bool Foam::ensightCase::options::separateCloud() const
{
    return separateCloud_;
}


void Foam::ensightCase::options::separateCloud(bool b)
{
    separateCloud_ = b;
}


// ************************************************************************* //

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

#include "foamVtkAppendRawFormatter.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::foamVtkAppendRawFormatter::name_     = "append";
const char* Foam::foamVtkAppendRawFormatter::encoding_ = "raw";


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::foamVtkAppendRawFormatter::write
(
    const char* s,
    std::streamsize n
)
{
    os().write(s, n);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkAppendRawFormatter::foamVtkAppendRawFormatter(std::ostream& os)
:
    foamVtkFormatter(os)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamVtkAppendRawFormatter::~foamVtkAppendRawFormatter()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const char* Foam::foamVtkAppendRawFormatter::name() const
{
    return name_;
}


const char* Foam::foamVtkAppendRawFormatter::encoding() const
{
    return encoding_;
}


void Foam::foamVtkAppendRawFormatter::writeSize(const uint64_t val)
{
    write(reinterpret_cast<const char*>(&val), sizeof(uint64_t));
}


void Foam::foamVtkAppendRawFormatter::write(const uint8_t val)
{
    write(reinterpret_cast<const char*>(&val), sizeof(uint8_t));
}


void Foam::foamVtkAppendRawFormatter::write(const label val)
{
    // std::cerr<<"label is:" << sizeof(val) << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(label));
}


void Foam::foamVtkAppendRawFormatter::write(const float val)
{
    // std::cerr<<"float is:" << sizeof(val) << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(float));
}


void Foam::foamVtkAppendRawFormatter::write(const double val)
{
    // std::cerr<<"write double as float:" << val << '\n';
    float copy(val);
    write(copy);
}


void Foam::foamVtkAppendRawFormatter::flush()
{}


// ************************************************************************* //

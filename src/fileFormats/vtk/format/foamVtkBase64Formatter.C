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

#include "foamVtkBase64Formatter.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::foamVtkBase64Formatter::name_     = "binary";
const char* Foam::foamVtkBase64Formatter::encoding_ = "base64";


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::foamVtkBase64Formatter::write
(
    const char* s,
    std::streamsize n
)
{
    base64Layer::write(s, n);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkBase64Formatter::foamVtkBase64Formatter(std::ostream& os)
:
    foamVtkFormatter(os),
    base64Layer(os)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamVtkBase64Formatter::~foamVtkBase64Formatter()
{
    flush();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const char* Foam::foamVtkBase64Formatter::name() const
{
    return name_;
}


const char* Foam::foamVtkBase64Formatter::encoding() const
{
    return encoding_;
}


void Foam::foamVtkBase64Formatter::writeSize(const uint64_t val)
{
    write(reinterpret_cast<const char*>(&val), sizeof(uint64_t));
}


void Foam::foamVtkBase64Formatter::write(const uint8_t val)
{
    base64Layer::add(val);
}


void Foam::foamVtkBase64Formatter::write(const label val)
{
    // std::cerr<<"label is:" << sizeof(val) << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(label));
}


void Foam::foamVtkBase64Formatter::write(const float val)
{
    // std::cerr<<"float is:" << sizeof(val) << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(float));
}


void Foam::foamVtkBase64Formatter::write(const double val)
{
    // std::cerr<<"write double as float:" << val << '\n';
    float copy(val);
    write(copy);
}


void Foam::foamVtkBase64Formatter::flush()
{
    if (base64Layer::close())
    {
        os().put('\n');
    }
}


// ************************************************************************* //

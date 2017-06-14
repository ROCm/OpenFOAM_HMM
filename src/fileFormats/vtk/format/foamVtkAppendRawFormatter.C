/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::vtk::appendRawFormatter::name_     = "append";
const char* Foam::vtk::appendRawFormatter::encoding_ = "raw";

const Foam::vtk::outputOptions
Foam::vtk::appendRawFormatter::opts_(formatType::APPEND_BINARY);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::vtk::appendRawFormatter::write
(
    const char* s,
    std::streamsize n
)
{
    os().write(s, n);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::appendRawFormatter::appendRawFormatter(std::ostream& os)
:
    formatter(os)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::appendRawFormatter::~appendRawFormatter()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::vtk::outputOptions&
Foam::vtk::appendRawFormatter::opts() const
{
    return opts_;
}


const char* Foam::vtk::appendRawFormatter::name() const
{
    return name_;
}


const char* Foam::vtk::appendRawFormatter::encoding() const
{
    return encoding_;
}


void Foam::vtk::appendRawFormatter::writeSize(const uint64_t nBytes)
{
    write(reinterpret_cast<const char*>(&nBytes), sizeof(uint64_t));
}


void Foam::vtk::appendRawFormatter::write(const uint8_t val)
{
    write(reinterpret_cast<const char*>(&val), sizeof(uint8_t));
}


void Foam::vtk::appendRawFormatter::write(const label val)
{
    // std::cerr<<"label:" << sizeof(val) << "=" << val << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(label));
}


void Foam::vtk::appendRawFormatter::write(const float val)
{
    // std::cerr<<"float:" << sizeof(val) << "=" << val << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(float));
}


void Foam::vtk::appendRawFormatter::write(const double val)
{
    // std::cerr<<"double as float=" << val << '\n';
    float copy(val);
    write(copy);
}


void Foam::vtk::appendRawFormatter::flush()
{/*nop*/}


// ************************************************************************* //

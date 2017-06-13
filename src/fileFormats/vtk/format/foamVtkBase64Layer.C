/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "foamVtkBase64Layer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::vtk::foamVtkBase64Layer::encoding_ = "base64";


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::vtk::foamVtkBase64Layer::write
(
    const char* s,
    std::streamsize n
)
{
    base64Layer::write(s, n);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::foamVtkBase64Layer::foamVtkBase64Layer(std::ostream& os)
:
    formatter(os),
    base64Layer(os)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::foamVtkBase64Layer::~foamVtkBase64Layer()
{
    base64Layer::close();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const char* Foam::vtk::foamVtkBase64Layer::encoding() const
{
    return encoding_;
}


void Foam::vtk::foamVtkBase64Layer::writeSize(const uint64_t nBytes)
{
    write(reinterpret_cast<const char*>(&nBytes), sizeof(uint64_t));
}


void Foam::vtk::foamVtkBase64Layer::write(const uint8_t val)
{
    base64Layer::add(val);
}


void Foam::vtk::foamVtkBase64Layer::write(const label val)
{
    // std::cerr<<"label:" << sizeof(val) << "=" << val << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(label));
}


void Foam::vtk::foamVtkBase64Layer::write(const float val)
{
    // std::cerr<<"float:" << sizeof(val) << "=" << val << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(float));
}


void Foam::vtk::foamVtkBase64Layer::write(const double val)
{
    // std::cerr<<"double as float=" << val << '\n';
    float copy(val);
    write(copy);
}


void Foam::vtk::foamVtkBase64Layer::flush()
{
    base64Layer::close();
}


std::size_t Foam::vtk::foamVtkBase64Layer::encodedLength
(
    std::size_t n
) const
{
    return base64Layer::encodedLength(n);
}


// ************************************************************************* //

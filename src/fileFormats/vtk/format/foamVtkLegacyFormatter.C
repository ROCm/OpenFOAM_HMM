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

#include "foamVtkLegacyFormatter.H"
#include "endian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::foamVtkLegacyFormatter::name_ = "BINARY";


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::foamVtkLegacyFormatter::write
(
    const char* s,
    std::streamsize n
)
{
    os().write(s, n);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkLegacyFormatter::foamVtkLegacyFormatter(std::ostream& os)
:
    foamVtkFormatter(os)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamVtkLegacyFormatter::~foamVtkLegacyFormatter()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const char* Foam::foamVtkLegacyFormatter::name() const
{
    return name_;
}


const char* Foam::foamVtkLegacyFormatter::encoding() const
{
    return name_;
}


void Foam::foamVtkLegacyFormatter::writeSize(const uint64_t)
{}


void Foam::foamVtkLegacyFormatter::write(const uint8_t val)
{
    // Legacy can only handle 32-bit integers.
    // Nonetheless promote to 'label' (32 or 64 bit) and deal with it later
    label copy(val);
    write(copy);
}


void Foam::foamVtkLegacyFormatter::write(const label val)
{
    // std::cerr<<"label is:" << sizeof(val) << '\n';

    // Not entirely correct: the legacy format only supports 32-bit integers.
    // Either limit size for 64-bit label, or simply do not support for 64-bit.
#ifdef WM_LITTLE_ENDIAN
# if WM_LABEL_SIZE == 32
    uint32_t swapped = endian::swap32(val);
    write(reinterpret_cast<const char*>(&swapped), sizeof(uint32_t));
# elif WM_LABEL_SIZE == 64
    uint64_t swapped = endian::swap64(val);
    write(reinterpret_cast<const char*>(&swapped), sizeof(uint64_t));
#endif
#else
    write(reinterpret_cast<const char*>(&val), sizeof(label));
#endif
}


void Foam::foamVtkLegacyFormatter::write(const float val)
{
    // std::cerr<<"float is:" << sizeof(val) << '\n';

#ifdef WM_LITTLE_ENDIAN
    // De-reference in two stages to avoid the warning
    //     dereferencing type-punned pointer will break strict-aliasing rules
    //     [-Wstrict-aliasing]

    const uint32_t* ptr = reinterpret_cast<const uint32_t*>(&val);
    uint32_t swapped = endian::swap32(*ptr);
    write(reinterpret_cast<const char*>(&swapped), sizeof(uint32_t));
#else
    write(reinterpret_cast<const char*>(&val), sizeof(float));
#endif
}


void Foam::foamVtkLegacyFormatter::write(const double val)
{
    // Legacy cannot support Float64 anyhow.
    // std::cerr<<"write double as float:" << val << '\n';
    float copy(val);
    write(copy);
}


void Foam::foamVtkLegacyFormatter::flush()
{
    os()<< '\n';
}


// ************************************************************************* //

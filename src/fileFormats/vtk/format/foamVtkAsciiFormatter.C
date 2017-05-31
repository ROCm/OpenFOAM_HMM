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

#include "foamVtkAsciiFormatter.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::foamVtkOutput::asciiFormatter::name_ = "ascii";

const Foam::foamVtkOutput::outputOptions
Foam::foamVtkOutput::asciiFormatter::opts_(formatType::INLINE_ASCII);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::foamVtkOutput::asciiFormatter::next()
{
    if (pos_ == 6)
    {
        os()<< '\n';
        pos_ = 0;
    }
    else if (pos_)
    {
        os()<< ' ';
    }
    ++pos_;
}


inline void Foam::foamVtkOutput::asciiFormatter::done()
{
    if (pos_)
    {
        os()<< '\n';
    }
    pos_ = 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkOutput::asciiFormatter::asciiFormatter(std::ostream& os)
:
    formatter(os),
    pos_(0)
{}


Foam::foamVtkOutput::asciiFormatter::asciiFormatter
(
    std::ostream& os,
    unsigned precision
)
:
    formatter(os),
    pos_(0)
{
    os.precision(precision);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamVtkOutput::asciiFormatter::~asciiFormatter()
{
    done();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::foamVtkOutput::outputOptions&
Foam::foamVtkOutput::asciiFormatter::opts() const
{
    return opts_;
}


const char* Foam::foamVtkOutput::asciiFormatter::name() const
{
    return name_;
}


const char* Foam::foamVtkOutput::asciiFormatter::encoding() const
{
    return name_;
}


void Foam::foamVtkOutput::asciiFormatter::writeSize(const uint64_t ignored)
{/*nop*/}


void Foam::foamVtkOutput::asciiFormatter::write(const uint8_t val)
{
    next();
    os()<< int(val);
}


void Foam::foamVtkOutput::asciiFormatter::write(const label val)
{
    next();
    os()<< val;
}


void Foam::foamVtkOutput::asciiFormatter::write(const float val)
{
    next();
    os()<< val;
}


void Foam::foamVtkOutput::asciiFormatter::write(const double val)
{
    next();
    os()<< float(val);
}


void Foam::foamVtkOutput::asciiFormatter::flush()
{
    done();
}


std::size_t
Foam::foamVtkOutput::asciiFormatter::encodedLength(std::size_t ignored) const
{
    return 0;
}


// ************************************************************************* //

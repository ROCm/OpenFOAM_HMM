/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "ensightFile.H"
#include "error.H"
#include "UList.H"
#include <cstring>
#include <sstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::ensightFile::allowUndef_ = false;

Foam::scalar Foam::ensightFile::undefValue_ = Foam::floatScalarVGREAT;

// Default is width 8
Foam::string Foam::ensightFile::mask_   = "********";
Foam::string Foam::ensightFile::dirFmt_ = "%08d";

const char* const Foam::ensightFile::coordinates = "coordinates";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::string Foam::ensightFile::mask()
{
    return mask_;
}


Foam::string Foam::ensightFile::subDir(const label n)
{
    char buf[32];

    sprintf(buf, dirFmt_.c_str(), n);
    return buf;
}


void Foam::ensightFile::subDirWidth(const label n)
{
    // enforce max limit to avoid buffer overflow in subDir()
    if (n < 1 || n > 31)
    {
        return;
    }

    // appropriate printf format
    std::ostringstream oss;
    oss << "%0" << n << "d";
    dirFmt_ = oss.str();

    // set mask accordingly
    mask_.resize(n, '*');
}


Foam::label Foam::ensightFile::subDirWidth()
{
    return mask_.size();
}


bool Foam::ensightFile::isUndef(const UList<scalar>& field)
{
    for (const scalar& val : field)
    {
        if (std::isnan(val))
        {
            return true;
        }
    }

    return true;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightFile::initialize()
{
    // ascii formatting specs
    setf
    (
        ios_base::scientific,
        ios_base::floatfield
    );
    precision(5);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightFile::ensightFile
(
    const fileName& pathname,
    IOstream::streamFormat format
)
:
    OFstream(ensight::FileName(pathname), format)
{
    initialize();
}


Foam::ensightFile::ensightFile
(
    const fileName& path,
    const fileName& name,
    IOstream::streamFormat format
)
:
    OFstream(path/ensight::FileName(name), format)
{
    initialize();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ensightFile::allowUndef()
{
    return allowUndef_;
}


bool Foam::ensightFile::allowUndef(bool enabled)
{
    bool old = allowUndef_;
    allowUndef_ = enabled;
    return old;
}


Foam::scalar Foam::ensightFile::undefValue(const scalar value)
{
    // enable its use too
    allowUndef_ = true;

    scalar old = undefValue_;
    undefValue_ = value;
    return old;
}


Foam::Ostream& Foam::ensightFile::writeString(const char* str)
{
    // Output 80 chars, but allocate for trailing nul character
    // to avoid -Wstringop-truncation warnings/errors.

    char buf[80+1];
    strncpy(buf, str, 80); // max 80 chars or padded with nul if smaller

    if (format() == IOstream::BINARY)
    {
        write(buf, 80);
    }
    else
    {
        buf[79] = 0;  // max 79 in ASCII, ensure it is indeed nul-terminated
        stdStream() << buf;
    }

    return *this;
}


Foam::Ostream& Foam::ensightFile::writeString(const std::string& str)
{
    return writeString(str.c_str());
}


Foam::Ostream& Foam::ensightFile::write(const char* str)
{
    return writeString(str);
}


Foam::Ostream& Foam::ensightFile::write(const word& str)
{
    return writeString(str);
}


Foam::Ostream& Foam::ensightFile::write(const string& str)
{
    return writeString(str);
}


Foam::Ostream& Foam::ensightFile::write
(
    const char* buf,
    std::streamsize count
)
{
    stdStream().write(buf, count);
    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const int32_t val)
{
    if (format() == IOstream::BINARY)
    {
        write
        (
            reinterpret_cast<const char *>(&val),
            sizeof(int32_t)
        );
    }
    else
    {
        stdStream().width(10);
        stdStream() << val;
    }

    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const int64_t val)
{
    int32_t ivalue(narrowInt32(val));

    return write(ivalue);
}


Foam::Ostream& Foam::ensightFile::write(const floatScalar val)
{
    if (format() == IOstream::BINARY)
    {
        write
        (
            reinterpret_cast<const char *>(&val),
            sizeof(floatScalar)
        );
    }
    else
    {
        stdStream().width(12);
        stdStream() << val;
    }

    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const doubleScalar val)
{
    float fvalue(narrowFloat(val));

    return write(fvalue);
}


Foam::Ostream& Foam::ensightFile::write
(
    const label value,
    const label fieldWidth
)
{
    if (format() == IOstream::BINARY)
    {
        write(value);
    }
    else
    {
        stdStream().width(fieldWidth);
        stdStream() << value;
    }

    return *this;
}


void Foam::ensightFile::newline()
{
    if (format() == IOstream::ASCII)
    {
        stdStream() << nl;
    }
}


Foam::Ostream& Foam::ensightFile::writeUndef()
{
    write(undefValue_);
    return *this;
}


Foam::Ostream& Foam::ensightFile::writeKeyword(const keyType& key)
{
    if (allowUndef_)
    {
        writeString(key + " undef");
        newline();
        write(undefValue_);
        newline();
    }
    else
    {
        writeString(key);
        newline();
    }

    return *this;
}


Foam::Ostream& Foam::ensightFile::writeBinaryHeader()
{
    if (format() == IOstream::BINARY)
    {
        writeString("C Binary");
    }

    return *this;
}


//
// Convenience Output Methods
//

void Foam::ensightFile::beginPart(const label index)
{
    writeString("part");
    newline();
    write(index+1); // Ensight starts with 1
    newline();
}


void Foam::ensightFile::beginParticleCoordinates(const label nparticles)
{
    writeString("particle coordinates");
    newline();
    write(nparticles, 8); // unusual width
    newline();
}


void Foam::ensightFile::writeLabels(const UList<label>& list)
{
    for (const label val : list)
    {
        write(val);
        newline();
    }
}


void Foam::ensightFile::writeList(const UList<label>& field)
{
    for (const label val : field)
    {
        write(scalar(val));
        newline();
    }
}


void Foam::ensightFile::writeList(const UList<scalar>& field)
{
    for (const scalar val : field)
    {
        if (std::isnan(val))
        {
            writeUndef();
        }
        else
        {
            write(val);
        }
        newline();
    }
}


// ************************************************************************* //

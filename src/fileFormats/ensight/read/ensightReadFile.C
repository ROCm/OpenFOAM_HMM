/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "ensightReadFile.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::IOstreamOption::streamFormat
Foam::ensightReadFile::detectBinaryHeader(const fileName& pathname)
{
    IOstreamOption::streamFormat fmt(IOstreamOption::BINARY);

    // Detect BINARY vs ASCII by testing for initial "(C|Fortran) Binary"
    {
        IFstream ifs(pathname, IOstreamOption::BINARY);

        if (!ifs.good())
        {
            FatalErrorInFunction
                << "Cannot read file " << ifs.name() << nl
                << exit(FatalError);
        }

        istream& iss = ifs.stdStream();

        // Binary string is *exactly* 80 characters
        string buf(size_t(80), '\0');
        iss.read(&buf[0], 80);

        if (!iss)
        {
            // Truncated?
            buf.erase(iss.gcount());
        }

        // Truncate at the first embedded '\0'
        const auto endp = buf.find('\0');
        if (endp != std::string::npos)
        {
            buf.erase(endp);
        }

        // ASCII if it does not contain "C Binary"
        if (!buf.contains("Binary") && !buf.contains("binary"))
        {
            fmt = IOstreamOption::ASCII;
        }
    }

    return fmt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightReadFile::ensightReadFile
(
    const fileName& pathname
)
:
    IFstream(pathname, ensightReadFile::detectBinaryHeader(pathname))
{}


Foam::ensightReadFile::ensightReadFile
(
    const fileName& pathname,
    IOstreamOption::streamFormat fmt
)
:
    IFstream(pathname, fmt)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Istream& Foam::ensightReadFile::read
(
    char* buf,
    std::streamsize count
)
{
    stdStream().read(buf, count);
    syncState();
    return *this;
}


Foam::Istream& Foam::ensightReadFile::read(string& value)
{
    if (format() == IOstreamOption::BINARY)
    {
        auto& iss = stdStream();

        // Binary string is *exactly* 80 characters
        value.resize(80, '\0');
        iss.read(&value[0], 80);

        syncState();

        if (!iss)
        {
            // Truncated - could also exit here, but no real advantage
            value.erase(iss.gcount());
        }

        // Truncate at the first embedded '\0'
        auto endp = value.find('\0');

        if (endp != std::string::npos)
        {
            value.erase(endp);
        }

        // May have been padded with trailing spaces - remove those
        endp = value.find_last_not_of(" \t\f\v\n\r");

        if (endp != std::string::npos)
        {
            value.erase(endp + 1);
        }
    }
    else
    {
        value.clear();
        while (value.empty() && !eof())
        {
            getLine(value);
        }
    }

    return *this;
}


Foam::Istream& Foam::ensightReadFile::read(label& value)
{
    int ivalue;

    if (format() == IOstreamOption::BINARY)
    {
        read
        (
            reinterpret_cast<char*>(&ivalue),
            sizeof(ivalue)
        );
    }
    else
    {
        stdStream() >> ivalue;
        syncState();
    }

    value = ivalue;
    return *this;
}


Foam::Istream& Foam::ensightReadFile::read(float& value)
{
    if (format() == IOstreamOption::BINARY)
    {
        read
        (
            reinterpret_cast<char*>(&value),
            sizeof(value)
        );
    }
    else
    {
        stdStream() >> value;
        syncState();
    }

    return *this;
}


Foam::Istream& Foam::ensightReadFile::read(double& value)
{
    float fvalue;

    if (format() == IOstreamOption::BINARY)
    {
        read
        (
            reinterpret_cast<char*>(&fvalue),
            sizeof(fvalue)
        );
    }
    else
    {
        stdStream() >> fvalue;
        syncState();
    }

    value = fvalue;
    return *this;
}


Foam::Istream& Foam::ensightReadFile::readKeyword(string& key)
{
    read(key);
    return *this;
}


Foam::Istream& Foam::ensightReadFile::readBinaryHeader()
{
    if (format() == IOstreamOption::BINARY)
    {
        string buffer;
        read(buffer);
    }

    return *this;
}


// ************************************************************************* //

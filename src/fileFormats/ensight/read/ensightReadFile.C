/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
    return *this;
}


Foam::Istream& Foam::ensightReadFile::read(string& value)
{
    if (format() == IOstream::BINARY)
    {
        auto& iss = stdStream();

        // Binary string is *exactly* 80 characters
        value.resize(80, '\0');
        iss.read(&value[0], 80);

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

    if (format() == IOstream::BINARY)
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
    }

    value = ivalue;
    return *this;
}


Foam::Istream& Foam::ensightReadFile::read(scalar& value)
{
    float fvalue;

    if (format() == IOstream::BINARY)
    {
        read
        (
            reinterpret_cast<char*>(&fvalue),
            sizeof(fvalue)
        );

        value = fvalue;
    }
    else
    {
        stdStream() >> value;
    }

    return *this;
}


Foam::Istream& Foam::ensightReadFile::readKeyword(string& key)
{
    read(key);
    return *this;
}


Foam::Istream& Foam::ensightReadFile::readBinaryHeader()
{
    if (format() == IOstream::BINARY)
    {
        string buffer;
        read(buffer);
    }

    return *this;
}


// ************************************************************************* //

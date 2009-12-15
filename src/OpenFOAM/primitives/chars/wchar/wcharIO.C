/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Reads a char from an input stream, for a given version
    number and File format. If an ascii File is being read, then the line
    numbers are counted and an erroneous read is reported.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "wchar.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const wchar_t wc)
{
    os.write(char(wc));

    os.check("Ostream& operator<<(Ostream&, const wchar_t)");
    return os;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const wchar_t* ws)
{
    if (ws)
    {
        for (const wchar_t* p = ws; *p; ++p)
        {
            os.write(char(*p));
        }
    }
    os.check("Ostream& operator<<(Ostream&, const wchar_t*)");
    return os;
}


// ************************************************************************* //

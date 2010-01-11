/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "wchar.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const wchar_t wc)
{
    if (!(wc & ~0x0000007F))
    {
        // 0x00000000 - 0x0000007F: (1-byte output)
        // 0xxxxxxx
        os.write(char(wc));
    }
    else if (!(wc & ~0x000007FF))
    {
        // 0x00000080 - 0x000007FF: (2-byte output)
        // 110bbbaa 10aaaaaa
        os.write(char(0xC0 | ((wc >> 6) & 0x1F)));
        os.write(char(0x80 | ((wc) & 0x3F)));
    }
    else if (!(wc & ~0x0000FFFF))
    {
        // 0x00000800 - 0x0000FFFF: (3-byte output)
        // 1110bbbb 10bbbbaa 10aaaaaa
        os.write(char(0xE0 | ((wc >> 12) & 0x0F)));
        os.write(char(0x80 | ((wc >> 6) & 0x3F)));
        os.write(char(0x80 | ((wc) & 0x3F)));
    }
    else if (!(wc & ~0x001FFFFF))
    {
        // 0x00010000 - 0x001FFFFF: (4-byte output)
        // 11110ccc 10ccbbbb 10bbbbaa 10aaaaaa
        os.write(char(0xF0 | ((wc >> 18) & 0x07)));
        os.write(char(0x80 | ((wc >> 12) & 0x3F)));
        os.write(char(0x80 | ((wc >> 6) & 0x3F)));
        os.write(char(0x80 | ((wc) & 0x3F)));
    }
    else if (!(wc & ~0x03FFFFFF))
    {
        // 0x00200000 - 0x03FFFFFF: (5-byte output)
        // 111110dd 10cccccc 10ccbbbb 10bbbbaa 10aaaaaa
        os.write(char(0xF8 | ((wc >> 24) & 0x03)));
        os.write(char(0x80 | ((wc >> 18) & 0x3F)));
        os.write(char(0x80 | ((wc >> 12) & 0x3F)));
        os.write(char(0x80 | ((wc >> 6) & 0x3F)));
        os.write(char(0x80 | ((wc) & 0x3F)));
    }
    else if (!(wc & ~0x7FFFFFFF))
    {
        // 0x04000000 - 0x7FFFFFFF: (6-byte output)
        // 1111110d 10dddddd 10cccccc 10ccbbbb 10bbbbaa 10aaaaaa
        os.write(char(0xFC | ((wc >> 30) & 0x01)));
        os.write(char(0x80 | ((wc >> 24) & 0x3F)));
        os.write(char(0x80 | ((wc >> 18) & 0x3F)));
        os.write(char(0x80 | ((wc >> 12) & 0x3F)));
        os.write(char(0x80 | ((wc >> 6) & 0x3F)));
        os.write(char(0x80 | ((wc) & 0x3F)));
    }
    else
    {
        // according to man page utf8(7)
        // the Unicode standard specifies no characters above 0x0010FFFF,
        // so Unicode characters can only be up to four bytes long in UTF-8.

        // report anything unknown/invalid as replacement character U+FFFD
        os.write(char(0xEF));
        os.write(char(0xBF));
        os.write(char(0xBD));
    }

    os.check("Ostream& operator<<(Ostream&, const wchar_t)");
    return os;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const wchar_t* ws)
{
    if (ws)
    {
        for (const wchar_t* p = ws; *p; ++p)
        {
            os  << *p;
        }
    }

    return os;
}


// ************************************************************************* //

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

Description

\*---------------------------------------------------------------------------*/

#include "UIBufStream.H"
#include "UOBufStream.H"
#include "wordList.H"
#include "IOstreams.H"
#include "argList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Buffer storage
    DynamicList<char> storage(1000);

    UOBufStream obuf(storage);
    obuf << 1002 << " " << "abcd" << " " << "def" << " " << 3.14159 << ";\n";

    Info<<"formatted: " << obuf.size() << " chars" << endl;

    // Match size
    storage.resize(obuf.size());

    Info<<"as string: " << string(storage.cdata(), storage.size()) << endl;

    // Attach input buffer - could also do without previous resize

    UIBufStream ibuf(storage, storage.size());

    token t;

    while (ibuf.good())
    {
        ibuf >> t;
        if (t.good())
        {
            Info<<"token: " << t << endl;
        }
    }

    Info<< nl << "Repeat..." << endl;
    ibuf.rewind();

    while (ibuf.good())
    {
        ibuf >> t;
        if (t.good())
        {
            Info<<"token: " << t << endl;
        }
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

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

#include "OListStream.H"
#include "wordList.H"
#include "IOstreams.H"
#include "argList.H"

using namespace Foam;

Ostream& toString(Ostream& os, const UList<char>& list)
{
    os << '"';
    for (const char c : list)
    {
        os << c;
    }
    os << '"';

    return os;
}


void printInfo(const OListStream& buf)
{
    Info<< nl << buf.size() << " chars (" << buf.capacity() << " capacity) ";
    toString(Info, buf.list()) << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Buffer storage
    DynamicList<char> storage(8);

    OListStream obuf(std::move(storage));
    obuf << 1002 << " " << "abcd" << " " << "def" << " " << 3.14159 << ";\n";

    printInfo(obuf);

    obuf.rewind();
    obuf << 100;

    printInfo(obuf);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "UListStream.H"
#include "wordList.H"
#include "IOstreams.H"
#include "argList.H"

#include <sstream>
#include <vector>

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


Ostream& toString(Ostream& os, const std::vector<char>& list)
{
    os << '"';
    for (const char c : list)
    {
        os << c;
    }
    os << '"';

    return os;
}


template<class BufType>
void printInfo(const BufType& buf)
{
    Info<< nl << "=========================" << endl;
    buf.print(Info);
    toString(Info, buf.list());
    Info<< nl << "=========================" << endl;
}


template<>
void printInfo(const UList<char>& buf)
{
    Info<< nl << "=========================" << endl;
    toString(Info, buf);
    Info<< nl << "=========================" << endl;
}


void printTokens(Istream& is)
{
    label count = 0;
    token t;
    while (is.good())
    {
        is >> t;
        if (t.good())
        {
            ++count;
            Info<<"token: " << t << endl;
        }
    }

    Info<< count << " tokens" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Buffer storage
    DynamicList<char> storage(1000);

    UOListStream obuf(storage);
    obuf << 1002 << "\n" << "abcd" << "\n" << "def" << "\n" << 3.14159 << ";\n";

    obuf.print(Info);

    // Match size
    storage.resize(obuf.size());

    Info<<"as string: " << string(storage.cdata(), storage.size()) << endl;

    // Attach input buffer - could also do without previous resize
    {
        UIListStream ibuf(storage);

        printTokens(ibuf);

        Info<< nl << "Repeat..." << endl;
        ibuf.rewind();

        printTokens(ibuf);
    }


    // Attach input buffer - could also do without previous resize
    {
        Info<< "parse as std::istream\n";

        uiliststream is(storage.cdata(), storage.size());

        string tok;

        while (std::getline(is, tok))
        {
            std::cerr << "tok: " << tok << nl;
        }

        Info<< nl << "Repeat..." << endl;

        is.rewind();
        while (std::getline(is, tok))
        {
            std::cerr << "tok: " << tok << nl;
        }
    }


    // Simple test using stl items only
    {
        std::vector<char> chars;
        chars.reserve(1024);

        for (int i=0; i < 10; ++i)
        {
            chars.push_back('A' + i);
            chars.push_back('0' + i);
            chars.push_back('\n');
        }

        Info<< "parse std::vector of char: ";
        toString(Info, chars);
        Info<< "----" << nl;

        uiliststream is(chars.data(), chars.size());
        string tok;
        std::cerr<< nl << "Parsed..." << nl;
        while (std::getline(is, tok))
        {
            std::cerr << "tok: " << tok << nl;
        }

        Info<< nl << "Repeat..." << endl;

        is.rewind();
        while (std::getline(is, tok))
        {
            std::cerr << "tok: " << tok << nl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

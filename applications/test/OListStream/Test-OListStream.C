/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "ListStream.H"
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


template<class BufType>
void printInfo(const BufType& buf)
{
    Info<< nl << "=========================" << endl;
    buf.print(Info);
    toString(Info, buf.list());
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


// Generate some dictionary-like content
template<class OS>
void outputDict(OS& os)
{
    os.beginBlock("testDict");
    os.writeEntry("bool",   "false");
    os.writeEntry("scalar", 3.14159);
    os.endBlock();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Buffer storage
    DynamicList<char> storage(16);

    OListStream obuf(std::move(storage));

    obuf.setBlockSize(100);

    printInfo(obuf);

    // Fill with some content
    for (label i = 0; i < 50; ++i)
    {
        obuf<< 1002 << " " << "abcd" << " "
            << "def" << " " << 3.14159 << ";\n";
    }

    printInfo(obuf);

    obuf.rewind();
    printInfo(obuf);

    for (label i=0; i < 10; ++i)
    {
        obuf << "item" << i << "\n";
    }

    printInfo(obuf);

    obuf.shrink();

    Info<< "after shrink" << nl;
    printInfo(obuf);

    // Add some more
    for (label i=10; i < 15; ++i)
    {
        obuf << "more" << i << nl;
    }

    Info<< "appended more" << nl;
    printInfo(obuf);

    // Overwrite at some position
    obuf.stdStream().rdbuf()->pubseekpos(0.60 * obuf.size());
    obuf << "<" << nl << "OVERWRITE" << nl;

    Info<<"after overwrite" << nl;
    printInfo(obuf);

    Info<< "transfer contents to a List or IListStream" << nl;

    IListStream ibuf;
    // Reclaim data storage from OListStream -> IListStream
    {
        List<char> data;
        obuf.swap(data);
        ibuf.swap(data);
    }

    Info<<"original:";
    printInfo(obuf);

    Info<<"new input:" << nl;
    printInfo(ibuf);

    printTokens(ibuf);

    // Create from other storage types

    DynamicList<char> written;
    Info<< nl;
    {
        Info<<"create std::move(List)" << endl;
        List<char> list(16, 'A');

        Info<<"input:";
        toString(Info, list) << endl;

        OListStream buf1(std::move(list));

        Info<<"orig:";
        toString(Info, list) << endl;
        printInfo(buf1);

        for (label i = 0; i < 26; ++i)
        {
            buf1 << char('A' +i);
        }
        for (label i = 0; i < 26; ++i)
        {
            buf1 << char('a' +i);
        }

        Info<<"orig:";
        toString(Info, list) << endl;

        printInfo(buf1);

        // Move back to written
        buf1.swap(written);

        printInfo(buf1);
    }
    Info<<"'captured' content ";
    toString(Info, written);

    Info<< nl
        << "content size=" << written.size()
        << " capacity=" << written.capacity() << nl;


    Info<< nl << "Test dictionary" << nl;
    {
        OListStream os1;

        outputDict(os1);

        Info<< "Regular" << nl;
        printInfo(os1);
    }

    {
        OListStream os2;
        os2.indentSize(0);

        outputDict(os2);

        Info<< "Compact" << nl;
        printInfo(os2);
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

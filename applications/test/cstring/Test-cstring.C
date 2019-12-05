/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Test some string functionality

\*---------------------------------------------------------------------------*/

#include "CStringList.H"
#include "DynamicList.H"
#include "IOstreams.H"
#include "fileNameList.H"
#include "stringOps.H"
#include "stringList.H"
#include "wordList.H"
#include "SubStrings.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int print(int argc, char *argv[])
{
    Info<< "argc=" << argc << endl;
    for (int i=0; i<argc; ++i)
    {
        Info<< "  argv[" << i << "] = \"" << argv[i] << "\"" << endl;
    }
    return argc;
}

int print(const CStringList& cstrLst)
{
    return print(cstrLst.size(), cstrLst.strings());
}


// Using nullptr termination
int print(char *argv[])
{
    if (!argv)
    {
        Info<< "argv=null" << endl;
        return 0;
    }

    int i=0;
    while (argv[i])
    {
        Info<< "  argv[" << i << "] = \"" << argv[i] << "\"" << endl;
        ++i;
    }

    return i;
}


// Main program:

int main(int argc, char *argv[])
{
    DynamicList<string> dynlst;
    dynlst.reserve(16);

    dynlst.append("string1 with content");
    dynlst.append("string2 other content");
    dynlst.append("string3 more");
    dynlst.append("string4 done");

    {
        CStringList inC(dynlst);

        Info<< "input: " << dynlst << nl;
        print(inC);

        Info<< "null-terminated string list" << nl;
        print(inC.strings());

        Info<< "sublist: starting at " << inC.size()/2 << nl;
        print(inC.strings(inC.size()/2));

        Info<< nl;
    }

    {
        string testInput
        (
            " A   test   input   line with various spacing  "
            "  and  text to be split on whitespace  "
        );

        Info<< testInput << nl;
        SubStrings<string> args = stringOps::splitSpace(testInput);
        Info<< "split into " << args.size() << " args" << nl;

        CStringList inC(args);
        print(inC);

        Info<< "sublist: starting at " << inC.size()/2 << nl;
        print(inC.strings(inC.size()/2));

        Info<< nl;
    }

    Info<<"command-line with " << CStringList::count(argv) << " items"<< endl;

    print(argc, argv);

    Info<< nl;
    {
        dynlst.clear();
        for (int i=0; i<argc; ++i)
        {
            dynlst.append(argv[i]);
        }

        Info<< "input: " << dynlst << endl;
        CStringList inC(dynlst);
        inC.reset(dynlst);

        print(inC);
        Info<< "length: " << inC.length() << endl;
        std::cout.write(inC.data(), inC.length());
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Test some string functionality

\*---------------------------------------------------------------------------*/

#include "CStringList.H"
#include "DynamicList.H"
#include "IOstreams.H"
#include "fileNameList.H"
#include "stringList.H"
#include "wordList.H"

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


// Main program:

int main(int argc, char *argv[])
{
    DynamicList<string> dynlst;
    dynlst.reserve(16);

    dynlst.append("string1 with content");
    dynlst.append("string2 other content");
    dynlst.append("string3 done");

    {
        CStringList inC(dynlst);

        Info<< "input: " << dynlst << endl;
        print(inC);
    }

    Info<<"command-line with " << CStringList::count(argv) << " items"<< endl;

    print(argc, argv);
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

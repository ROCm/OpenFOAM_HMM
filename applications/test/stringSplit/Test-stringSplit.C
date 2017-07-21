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

Application
    Test-stringSplit

Description
    Test string splitting

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "stringOps.H"

using namespace Foam;

template<class String>
void printSplitting(const String& str, const char delimiter)
{
    auto split = stringOps::split(str, delimiter);

    Info<< "string {" << str.size() << " chars} = " << str << nl
        << split.size() << " elements {" << split.length() << " chars}"
        << nl;

    unsigned i = 0;
    for (const auto s : split)
    {
        Info<< "[" << i++ << "] {" << s.length() << " chars} = "
            << s.str() << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();

    argList args(argc, argv, false, true);

    if (args.size() <= 1 && args.options().empty())
    {
        args.printUsage();
    }

    for (label argi=1; argi < args.size(); ++argi)
    {
        printSplitting(args[argi], '/');
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //

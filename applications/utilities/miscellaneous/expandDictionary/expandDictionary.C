/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    expandDictionary

Description
    Read the dictionary provided as an argument, expand the macros etc. and
    write the resulting dictionary to standard output.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "IOobject.H"
#include "dictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Read the specified dictionary file, expand the macros etc. and write\n"
        "the resulting dictionary to standard output."
    );

    argList::noBanner();
    argList::noParallel();
    argList::validArgs.append("inputDict");
    argList args(argc, argv);

    const string dictName = args[1];

    IOobject::writeBanner(Info)
        <<"//\n// " << dictName << "\n//\n";

    dictionary(IFstream(dictName)(), true).write(Info, false);

    IOobject::writeDivider(Info);

    return 0;
}


// ************************************************************************* //

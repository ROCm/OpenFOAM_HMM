/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Application
    fileNameCleanTest

Description


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "SubList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "OSspecific.H"


using namespace Foam;

//
// * remove repeated slashes
//       /abc////def        -->   /abc/def
//
// * remove '/./'
//       /abc/def/./ghi/.   -->   /abc/def/./ghi
//       abc/def/./         -->   abc/def
//
// * remove '/../'
//       /abc/def/../ghi/jkl/nmo/..   -->   /abc/ghi/jkl
//       abc/../def/ghi/../jkl        -->   abc/../def/jkl
//
// * remove trailing '/'
//
bool fileNameCleanThis(fileName& This)
{
    // the top slash - we are never allowed to descend below this one
    register string::size_type top = This.find('/');

    // no slashes - nothing to do
    if (top == string::npos)
    {
        return false;
    }

    // start with the '/' found:
    register char prev = '/';
    register string::size_type nChar = top+1;
    register string::size_type maxLen = This.size();

    for
    (
        register string::size_type src = nChar;
        src < maxLen;
        /*nil*/
    )
    {
        register char c = This[src++];

        if (prev == '/')
        {
            // repeated '/' - skip it
            if (c == '/')
            {
                continue;
            }

            // could be '/./' or '/../'
            if (c == '.')
            {
                // found trailing '/.' - skip it
                if (src >= maxLen)
                {
                    continue;
                }


                // peek at the next character
                register char c1 = This[src];

                // found '/./' - skip it
                if (c1 == '/')
                {
                    src++;
                    continue;
                }

                // it is '/..' or '/../'
                if (c1 == '.' && (src+1 >= maxLen || This[src+1] == '/'))
                {
                    // find a candidate for the parent directory
                    // minimum of 3 characters:  '/x/../'
                    string::size_type parent =
                    (
                        nChar > 2
                      ? This.rfind('/', nChar-2)
                      : string::npos
                    );


                    // strip parent directory, provided that it isn't below
                    // the current top edit point
                    if (parent != string::npos && parent >= top)
                    {
                        nChar = parent + 1;
                        src += 2;
                        continue;
                    }

                    // bad resolution, eg 'abc/../../'
                    // retain it, but move the top to avoid it being
                    // considered a valid parent later
                    top += 3;
                }
            }
        }
        This[nChar++] = prev = c;
    }

    // remove trailing slash
    if (nChar > 1 && This[nChar-1] == '/')
    {
        nChar--;
    }

    This.resize(nChar);

    return (nChar != maxLen);
}


fileName fileNameClean(fileName& This)
{
    fileName fName(This);
    fileNameCleanThis(fName);

    return fName;
}


void printCleaning(fileName& pathName)
{
    Info<< "fileName = " << pathName << nl
        << "  path() = " << pathName.path() << nl
        << "  name() = " << pathName.name() << nl << nl;

    fileNameCleanThis(pathName);

    Info<< "cleaned  = " << pathName << nl
        << "  path() = " << pathName.path() << nl
        << "  name() = " << pathName.name() << nl << nl;

    IOobject::writeDivider(Info);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::validArgs.insert("fileName .. fileNameN");

    argList args(argc, argv, false, true);

    if (args.additionalArgs().empty())
    {
        args.printUsage();
    }

    if (args.optionFound("case"))
    {
        fileName pathName = args.option("case");
        Info<< "-case\n";

        printCleaning(pathName);
    }

    forAll(args.additionalArgs(), argI)
    {
        fileName pathName = args.additionalArgs()[argI];
        printCleaning(pathName);
    }


    Info<< "\nEnd" << endl;

    return 0;
}


// ************************************************************************* //

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

Description
    Test the tokenizing of various things
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "IStringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.insert("string .. stringN");
    argList::validOptions.insert("file", "name");

    argList args(argc, argv, false, true);

    forAll(args.additionalArgs(), argI)
    {
        const string& rawArg = args.additionalArgs()[argI];
        Info<< "input string: " << rawArg << nl;

        IStringStream is(rawArg);
        
        while (is.good())
        {
            token tok(is);
            Info<< "token: " << tok.info() << endl;
        }
        
        Info<< nl;
        IOobject::writeDivider(Info);
    }  
    
    
    if (args.optionFound("file"))
    {
        IFstream is(args.option("file"));
        
        Info<< "tokenizing file: " << args.option("file") << nl;

        while (is.good())
        {
            token tok(is);
            Info<< "token: " << tok.info() << endl;
        }
        
        Info<< nl;
        IOobject::writeDivider(Info);
    }

    return 0;
}

// ************************************************************************* //

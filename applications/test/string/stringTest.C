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

Description

\*---------------------------------------------------------------------------*/

#include "string.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    string test("$HOME kjhkjhkjh \" \\$HOME/tyetyery ${FOAM_RUN} \n ; hkjh ;$");

    Info<< test << endl;

    // test sub-strings via iterators
    string::const_iterator iter = test.end();
    string::const_iterator iter2 = test.end();
    string::size_type fnd = test.find('\\');

    if (fnd != string::npos)
    {
        iter  = test.begin() + fnd;
        iter2 = iter + 6;
    }

    Info<< "sub-string via iterators : >";
    while (iter != iter2)
    {
        Info<< *iter;
        iter++;
    }
    Info<< "<\n";

    Info<< string(test).replace("kj", "zzz") << endl;
    Info<< string(test).replace("kj", "") << endl;
    Info<< string(test).replaceAll("kj", "zzz") << endl;
    Info<< string(test).replaceAll("kj", "z") << endl;

    Info<< string(test).expand() << endl;

    string test2("~OpenFOAM/controlDict");
    Info<< test2 << " => " << test2.expand() << endl;

    string s;
    Sin.getLine(s);

    string s2(s.expand());

    cout<< "output string with " << s2.length() << " characters\n";
    cout<< "ostream<<  >" << s2 << "<\n";
    Info<< "Ostream<<  >" << s2 << "<\n";

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

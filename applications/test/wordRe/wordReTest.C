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

#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "List.H"
#include "Tuple2.H"
#include "wordRe.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    wordRe wre;
    std::string s1("this .* file");
    Foam::string s2("this .* file");
    const char * s3 = "this .* file";

    Info<< wordRe(s1).info() << endl;
    Info<< wordRe(s2, false).info() << endl;
    Info<< wordRe(s2).info() << endl;
    Info<< wordRe(s3, true).info() << endl;

    wre = "this .* file";
    Info<< wre.info() << endl;
    wre = s1;
    Info<< wre.info() << endl;
    wre.uncompile();
    Info<< wre.info() << " uncompiled" << endl;

    
    wre = "something";
    Info<< wre.info() << " before" << endl;
    wre.uncompile();
    Info<< wre.info() << " uncompiled" << endl;
    wre.compile(true);
    Info<< wre.info() << " after auto-detect" << endl;

    wre = "something .* value";
    Info<< wre.info() << " before" << endl;
    wre.uncompile();
    Info<< wre.info() << " uncompiled" << endl;
    wre.compile(true);
    Info<< wre.info() << " after auto-detect" << endl;
    wre.uncompile();
    Info<< wre.info() << " uncompiled" << endl;
    wre.recompile();
    Info<< wre.info() << " recompiled" << endl;

    IOobject::writeDivider(Info);

    List<Tuple2<wordRe, string> > rawList(IFstream("testRegexps")());
    Info<< "input list:" << rawList << endl;
    IOobject::writeDivider(Info);
    Info<< endl;

    forAll(rawList, elemI)
    {
        const wordRe& wre = rawList[elemI].first();
        const string& str = rawList[elemI].second();

        Info<< wre.info()
            << " equals:" << (wre == str)
            << "(" << wre.match(str, true) << ")"
            << " match:" << wre.match(str)
            << "  str=" << str
            << endl;
    }

    Info<< endl;

    return 0;
}


// ************************************************************************* //

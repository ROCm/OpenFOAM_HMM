/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "List.H"
#include "Tuple2.H"
#include "keyType.H"
#include "wordRe.H"
#include "wordRes.H"
#include "predicates.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    wordRe wre;
    std::string s1("this .* file");
    Foam::string s2("this .* file");
    const char * s3 = "this .* file";

    keyType keyre("x.*", true);

    wordReList wordrelist
    {
        {"this", wordRe::LITERAL},
        {"x.*", wordRe::REGEX},
        {"file[a-b]", wordRe::REGEX},
    };

    if (true)
    {
        Info<<"keyType: " << keyre << endl;

        keyType key2(std::move(keyre));

        Info<<"move construct: <" << keyre << "> <" << key2 << ">" << endl;

        keyre = std::move(key2);

        Info<<"move assign: <" << keyre << "> <" << key2 << ">" << endl;

        keyType key3;

        keyre.swap(key3);

        Info<<"swap: <" << keyre << "> <" << key3 << ">" << endl;

        keyre = std::move(key3);
        Info<<"move assign: <" << keyre << "> <" << key3 << ">" << endl;

        return 0;
    }

    if (false)
    {
        wordRe keyre("y.*", wordRe::REGEX);

        Info<<"wordRe: " << keyre << endl;

        wordRe key2(std::move(keyre));

        Info<<"keyTypes: " << keyre << " " << key2 << endl;

        keyre = std::move(key2);

        Info<<"keyTypes: " << keyre << " " << key2 << endl;

        wordRe key3;

        keyre.swap(key3);

        Info<<"keyTypes: <" << keyre << "> <" << key3 << ">" << endl;

        keyre = std::move(key3);
        Info<<"keyTypes: <" << keyre << "> <" << key3 << ">" << endl;

        return 0;
    }

    wordRes wrelist(wordrelist);

    Info<< "re-list:" << wrelist() << endl;
    Info<< "match this: " << wrelist("this") << endl;
    Info<< "match xyz: "  << wrelist("xyz") << endl;
    Info<< "match zyx: "  << wrelist("zyx") << endl;
    Info<< "match xyz: "  << wrelist.match("xyz") << endl;
    Info<< "match any: "  << predicates::always()("any junk") << endl;
    Info<< "keyre match: "  << keyre("xyz") << endl;
    Info<< "string match: "  << string("this").match("xyz") << endl;
    Info<< "string match: "  << string("x.*")("xyz") << endl;
    Info<< "string match: "  << string("x.*")(keyre) << endl;

    wordRe(s1, wordRe::DETECT).info(Info) << endl;
    wordRe(s2).info(Info) << endl;
    wordRe(s2, wordRe::DETECT).info(Info) << endl;
    wordRe(s3, wordRe::REGEX).info(Info) << endl;

    wre = "this .* file";

    Info<<"substring: " << wre.substr(4) << endl;

    wre.info(Info) << endl;
    wre = s1;
    wre.info(Info) << endl;
    wre.uncompile();
    wre.info(Info) << endl;

    wre = "something";
    wre.info(Info) << " before" << endl;
    wre.uncompile();
    wre.info(Info) << " uncompiled" << endl;
    wre.compile(wordRe::DETECT);
    wre.info(Info) << " after DETECT" << endl;
    wre.compile(wordRe::NOCASE);
    wre.info(Info) << " after NOCASE" << endl;
    wre.compile(wordRe::DETECT_NOCASE);
    wre.info(Info) << " after DETECT_NOCASE" << endl;

    wre = "something .* value";
    wre.info(Info) << " before" << endl;
    wre.uncompile();
    wre.info(Info) << " uncompiled" << endl;
    wre.compile(wordRe::DETECT);
    wre.info(Info) << " after DETECT" << endl;
    wre.uncompile();
    wre.info(Info) << " uncompiled" << endl;
    wre.compile();
    wre.info(Info) << " re-compiled" << endl;

    wre.set("something .* value", wordRe::LITERAL);
    wre.info(Info) << " set as LITERAL" << endl;

    IOobject::writeDivider(Info);

    List<Tuple2<wordRe, string>> rawList(IFstream("testRegexps")());
    Info<< "input list:" << rawList << endl;
    IOobject::writeDivider(Info) << endl;

    forAll(rawList, elemI)
    {
        const wordRe& wre = rawList[elemI].first();
        const string& str = rawList[elemI].second();

        wre.info(Info)
            << " equals:" << (wre == str)
            << "(" << wre.match(str, true) << ")"
            << " match:" << wre.match(str)
            << "  str=" << str
            << endl;

        wordRe wre2;
        wre2.set(wre, wordRe::NOCASE);

        wre2.info(Info)
            << " match:" << wre2.match(str)
            << "  str=" << str
            << endl;

    }

    Info<< endl;

    return 0;
}


// ************************************************************************* //

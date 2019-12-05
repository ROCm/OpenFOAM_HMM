/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Test word/regex

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "List.H"
#include "Tuple2.H"
#include "keyType.H"
#include "wordRes.H"
#include "predicates.H"
#include "Random.H"

using namespace Foam;


word typeOf(wordRe::compOption retval)
{
    if (wordRe::LITERAL == retval) return "(literal)";
    if (wordRe::UNKNOWN == retval) return "(unknown)";
    if (wordRe::REGEX == retval)   return "(regex)";
    return "";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    wordRe wre;
    std::string s1("this .* file");
    Foam::string s2("this .* file");
    const char * s3 = "this .* file";

    keyType keyre("x.*", keyType::REGEX);

    wordReList wordrelist
    {
        {"this", wordRe::LITERAL},
        {"x.*", wordRe::REGEX},
        {"file[a-b]", wordRe::REGEX},
        {"xvalues", wordRe::LITERAL},
        {"xv.*", wordRe::REGEX},
    };

    if (false)
    {
        Info<<"keyType: " << keyre << nl;

        keyType key2(std::move(keyre));

        Info<<"move construct: <" << keyre << "> <" << key2 << ">" << nl;

        keyre = std::move(key2);

        Info<<"move assign: <" << keyre << "> <" << key2 << ">" << nl;

        keyType key3;

        keyre.swap(key3);

        Info<<"swap: <" << keyre << "> <" << key3 << ">" << nl;

        keyre = std::move(key3);
        Info<<"move assign: <" << keyre << "> <" << key3 << ">" << nl;

        return 0;
    }

    if (false)
    {
        wordRe keyre("y.*", wordRe::REGEX);

        Info<<"wordRe: " << keyre << nl;

        wordRe key2(std::move(keyre));

        Info<<"keyTypes: " << keyre << " " << key2 << nl;

        keyre = std::move(key2);

        Info<<"keyTypes: " << keyre << " " << key2 << nl;

        wordRe key3;

        keyre.swap(key3);

        Info<<"keyTypes: <" << keyre << "> <" << key3 << ">" << nl;

        keyre = std::move(key3);
        Info<<"keyTypes: <" << keyre << "> <" << key3 << ">" << nl;

        return 0;
    }

    wordRes wres1(wordrelist);

    Info<< "re-list:" << wres1 << nl;
    Info<< "match this: " << wres1("this") << nl;
    Info<< "match xyz: "  << wres1("xyz") << nl;
    Info<< "match zyx: "  << wres1("zyx") << nl;
    Info<< "match xyz: "  << wres1.match("xyz") << nl;
    Info<< "match any: "  << predicates::always()("any junk") << nl;

    Info<< "match xvalues: "  << wres1.match("xvalues") << nl;
    Info<< "matched xvalues = "  << typeOf(wres1.matched("xvalues")) << nl;
    Info<< "matched xval = "  << typeOf(wres1.matched("xval")) << nl;
    Info<< "matched zyx = "  << typeOf(wres1.matched("zyx")) << nl;

    Info<< "keyre match: "  << keyre("xyz") << nl;
    Info<< "string match: "  << string("this").match("xyz") << nl;
    Info<< "string match: "  << string("x.*")("xyz") << nl;
    Info<< "string match: "  << string("x.*")(keyre) << nl;


    // Test uniq
    {
        Random rnd;
        const label last = wres1.size()-1;

        for (label i = 0; i < 8; ++i)
        {
            // Make a copy
            wordRe wre(wres1[rnd.position<label>(0,last)]);

            // Append
            wres1.append(wre);
        }

        // Add some entropy
        Foam::shuffle(wres1);

        Info<< nl
            << "Test uniq on " << wres1
            << "  ==  " << wordRes::uniq(wres1) << nl;

        // Inplace
        wres1.uniq();
        Info<< nl << "Inplace: " << wres1 << nl;
    }
    Info<< nl;

    wordRe(s1, wordRe::DETECT).info(Info) << nl;
    wordRe(s2).info(Info) << nl;
    wordRe(s2, wordRe::DETECT).info(Info) << nl;
    wordRe(s3, wordRe::REGEX).info(Info) << nl;

    wre = "this .* file";

    Info<<"substring: " << wre.substr(4) << nl;

    wre.info(Info) << nl;
    wre = s1;
    wre.info(Info) << nl;
    wre.uncompile();
    wre.info(Info) << nl;

    wre = "something";
    wre.info(Info) << " before" << nl;
    wre.uncompile();
    wre.info(Info) << " uncompiled" << nl;
    wre.compile(wordRe::DETECT);
    wre.info(Info) << " after DETECT" << nl;
    wre.compile(wordRe::ICASE);
    wre.info(Info) << " after ICASE" << nl;
    wre.compile(wordRe::DETECT_ICASE);
    wre.info(Info) << " after DETECT_ICASE" << nl;

    wre = "something .* value";
    wre.info(Info) << " before" << nl;
    wre.uncompile();
    wre.info(Info) << " uncompiled" << nl;
    wre.compile(wordRe::DETECT);
    wre.info(Info) << " after DETECT" << nl;
    wre.uncompile();
    wre.info(Info) << " uncompiled" << nl;
    wre.compile();
    wre.info(Info) << " re-compiled" << nl;

    wre.set("something .* value", wordRe::LITERAL);
    wre.info(Info) << " set as LITERAL" << nl;

    IOobject::writeDivider(Info);

    List<Tuple2<wordRe, string>> rawList(IFstream("testRegexps")());
    Info<< "input list:" << rawList << nl;
    IOobject::writeDivider(Info) << nl;

    forAll(rawList, elemI)
    {
        const wordRe& wre = rawList[elemI].first();
        const string& str = rawList[elemI].second();

        wre.info(Info)
            << " equals:" << (wre == str)
            << "(" << wre.match(str, true) << ")"
            << " match:" << wre.match(str)
            << "  str=" << str
            << nl;

        wordRe wre2;
        wre2.set(wre, wordRe::ICASE);

        wre2.info(Info)
            << " match:" << wre2.match(str)
            << "  str=" << str
            << nl;

    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

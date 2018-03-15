/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Tests for regular expressions

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "regExp.H"
#include "SubStrings.H"
#include "Switch.H"

using namespace Foam;


//- Simple structure for holding regex tests
struct regexTest
{
    bool expected;
    string pattern;
    string text;

    friend Istream& operator>>(Istream& is, struct regexTest& obj)
    {
        is.readBegin("struct");
        is >> obj.expected >> obj.pattern >> obj.text;
        is.readEnd("struct");
        is.check(FUNCTION_NAME);
        return is;
    }

    friend Ostream& operator<<(Ostream& os, const struct regexTest& obj)
    {
        os  << "( "
            << obj.expected << " " << obj.pattern << " " << obj.text
            << " )";
        return os;
    }
};

// Needed for list output. Just treat everything as unequal.
bool operator!=(const struct regexTest&, const struct regexTest&)
{
    return true;
}

// Simple output of match groups
static Ostream& operator<<(Ostream& os, const regExp::results_type& sm)
{
    for (std::smatch::size_type i = 1; i < sm.size(); ++i)
    {
        os << " " << sm.str(i);
    }

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    List<regexTest> rawList(IFstream("testRegexps")());
    Info<< "Test expressions:" << rawList << endl;
    IOobject::writeDivider(Info) << endl;

    regExp::results_type match;

    // Expect some failures:
    const bool throwingError = FatalError.throwExceptions();

    // Report matches:
    for (const auto& testseq : rawList)
    {
        const bool expected = testseq.expected;
        const string& pat = testseq.pattern;
        const string& str = testseq.text;

        Info<< "Test " << Switch(expected) << ": "
            << str << " =~ m/" << pat.c_str() << "/ == ";

        regExp re;

        try
        {
            re = pat;

            if (re.match(str, match))
            {
                Info<< "true";
                if (re.ngroups())
                {
                    Info<< " (" << re.ngroups() << " groups):" << match;
                }
            }
            else if (re.search(str))
            {
                Info<< "partial match";
            }
            else
            {
                Info<< "false";
            }
            Info<< endl;
        }
        catch (Foam::error& err)
        {
            Info<< "Caught FatalError " << err << nl << endl;
            continue;
        }

        if (false)
        {
            regExp re2(std::move(re));
            Info<<"move construct: " << re.exists() << "/" << re2.exists()
                << endl;

            re = std::move(re2);
            Info<<"move assign: " << re.exists() << "/" << re2.exists()
                << endl;

            re.swap(re2);
            Info<<"swap: " << re.exists() << "/" << re2.exists()
                << endl;
        }
    }

    Info<< nl << "test regExp(const char*) ..." << endl;
    string me("Mark");

    try
    {
        // Handling of null strings
        if (regExp(nullptr).match(me))
        {
            Info<< "fail - matched: " << me << endl;
        }
        else
        {
            Info<< "pass - null pointer is no expression" << endl;
        }
    }
    catch (Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    try
    {
        // Normal match
        if (regExp("[Mm]ar[ck]").match(me))
        {
            Info<< "pass - matched: " << me << endl;
        }
        else
        {
            Info<< "no match" << endl;
        }
    }
    catch (Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    try
    {
        // Match ignore case
        if (regExp("mar[ck]", true).match(me))
        {
            Info<< "pass - matched: " << me << endl;
        }
        else
        {
            Info<< "no match" << endl;
        }
    }
    catch (Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    try
    {
        // Embedded prefix for match ignore case
        if (regExp("(?i)mar[ck]").match(me))
        {
            Info<< "pass - matched: " << me << endl;
        }
        else
        {
            Info<< "no match" << endl;
        }
    }
    catch (Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    try
    {
        // Handling of empty expression
        if (regExp("").match(me))
        {
            Info<< "fail - matched: " << me << endl;
        }
        else
        {
            Info<< "pass - no match on empty expression" << endl;
        }
    }
    catch (Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    try
    {
        // Embedded prefix - but expression is empty
        if (regExp("(?i)").match(me))
        {
            Info<< "fail - matched: " << me << endl;
        }
        else
        {
            Info<< "pass - no match on empty expression" << endl;
        }
    }
    catch (Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    FatalError.throwExceptions(throwingError);

    Info<< "\nDone" << nl << endl;

    return 0;
}


// ************************************************************************* //

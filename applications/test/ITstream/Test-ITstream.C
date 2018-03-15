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

Description

\*---------------------------------------------------------------------------*/

#include "ListStream.H"
#include "UListStream.H"
#include "wordList.H"
#include "IOstreams.H"
#include "argList.H"
#include "ITstream.H"

using namespace Foam;

template<class T>
Ostream& toString(Ostream& os, const T& str)
{
    os << str;
    return os;
}


template<>
Ostream& toString(Ostream& os, const UList<char>& list)
{
    for (const char c : list)
    {
        os << c;
    }

    return os;
}


template<>
Ostream& toString(Ostream& os, const List<char>& list)
{
    for (const char c : list)
    {
        os << c;
    }

    return os;
}


void printTokens(Istream& is)
{
    label count = 0;
    token t;
    while (is.good())
    {
        is >> t;
        if (t.good())
        {
            ++count;
            Info<<"token: " << t << endl;
        }
    }

    Info<< count << " tokens" << endl;
}


template<class BUF>
void doTest(const string& name, const BUF& input, bool verbose=false)
{
    Info<<"test " << name.c_str() << ":" << nl
        <<"====" << nl;
    toString(Info, input)
        << nl
        <<"====" << nl << endl;

    ITstream its(name, input);
    Info<< "got " << its.size() << " tokens - index at "
        << its.tokenIndex() << endl;

    if (verbose)
    {
        for (const token& tok : its)
        {
            Info<< "    " << tok.info() << nl;
        }
        Info<< nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    const char* charInput =
        "( const char input \"string\" to tokenize )"
        "List<label> 5(0 1 2 3 4);";

    string stringInput("( string     ;    input \"string\" to tokenize )");

    List<char> listInput(stringInput.cbegin(), stringInput.cend());

    doTest("char*", charInput, true);
    doTest("string", stringInput, true);
    doTest("List<char>", listInput, true);

    reverse(listInput);
    doTest("List<char>", listInput, true);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

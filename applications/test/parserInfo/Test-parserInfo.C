/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Output some (expressions) parser information

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "List.H"
#include "fieldExprParser.H"
#include "patchExprParser.H"
#include "volumeExprParser.H"

using namespace Foam;

template<class Parser>
void printInformation
(
    Ostream& os,
    const word& name,
    const bool printNames,
    const bool printRules
)
{
    if (printNames)
    {
        os << nl << name << " tokenNames:" << nl;
        Parser::printTokenNames(os);
    }

    if (printRules)
    {
        os << nl << name << " rules:" << nl;
        Parser::printRules(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::addNote
    (
        "Display token names or rules for specified expression parser(s)."
        " Without options, displays everything."
    );

    argList::addBoolOption("rules", "Print parser rules");
    argList::addBoolOption("tokens", "Print token names");

    argList::addBoolOption("field", "Field expression parser");
    argList::addBoolOption("patch", "Patch expression parser");
    argList::addBoolOption("volume", "Volume expression parser");

    argList args(argc, argv);

    // Defaults
    const bool all = !args.count({"field", "patch", "volume"});
    const bool both = !args.count({"tokens", "rules"});

    const bool printNames = both || args.found("tokens");
    const bool printRules = both || args.found("rules");


    if (all || args.found("field"))
    {
        printInformation<Foam::expressions::fieldExpr::parser>
        (
            Info,
            "field",
            printNames,
            printRules
        );
    }

    if (all || args.found("patch"))
    {
        printInformation<Foam::expressions::patchExpr::parser>
        (
            Info,
            "patch",
            printNames,
            printRules
        );
    }

    if (all || args.found("volume"))
    {
        printInformation<Foam::expressions::volumeExpr::parser>
        (
            Info,
            "volume",
            printNames,
            printRules
        );
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

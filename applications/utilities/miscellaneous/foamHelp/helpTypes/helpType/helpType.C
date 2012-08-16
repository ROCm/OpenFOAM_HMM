/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "helpType.H"
#include "doxygenXmlParser.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(helpType, 0);
    defineRunTimeSelectionTable(helpType, dictionary);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::helpType::displayDocOptions
(
    const string& searchStr,
    const bool exactMatch
) const
{
    const dictionary& docDict = debug::controlDict().subDict("Documentation");
    List<fileName> docDirs(docDict.lookup("doxyDocDirs"));

    label i = -1;
    forAll(docDirs, dirI)
    {
        if (isDir(docDirs[dirI].expand()))
        {
            i = dirI;
            break;
        }
    }

    if (i != -1)
    {
        Info<< "Found doxygen help in " << docDirs[i].c_str() << nl << endl;

        string docBrowser = getEnv("FOAM_DOC_BROWSER");
        if (docBrowser.empty())
        {
            docDict.lookup("docBrowser") >> docBrowser;
        }

        doxygenXmlParser parser
        (
            docDirs[i]/"../DTAGS",
            "tagfile",
            searchStr,
            exactMatch
        );

        Info<< "Valid types include:" << nl << SortableList<word>(parser.toc());
    }
    else
    {
        Info<< "No Doxygen sources found under search paths: "
            << docDirs << endl;
    }
}


void Foam::helpType::displayDoc
(
    const word& className,
    const string& searchStr,
    const bool exactMatch
) const
{
    const dictionary& docDict = debug::controlDict().subDict("Documentation");
    List<fileName> docDirs(docDict.lookup("doxyDocDirs"));

    label i = -1;
    forAll(docDirs, dirI)
    {
        if (isDir(docDirs[dirI].expand()))
        {
            i = dirI;
            break;
        }
    }

    if (i != -1)
    {
        Info<< "Found doxygen help in " << docDirs[i].c_str() << nl << endl;

        string docBrowser = getEnv("FOAM_DOC_BROWSER");
        if (docBrowser.empty())
        {
            docDict.lookup("docBrowser") >> docBrowser;
        }

        doxygenXmlParser parser
        (
            docDirs[i]/"../DTAGS",
            "tagfile",
            searchStr,
            exactMatch
        );

        if (debug)
        {
            Info<< parser;
        }

        if (parser.found(className))
        {
            fileName docFile
            (
                docDirs[i]/parser.subDict(className).lookup("filename")
            );

            // can use FOAM_DOC_BROWSER='application file://%f' if required
            docBrowser.replaceAll("%f", docFile);

            fileName classFolder(parser.subDict(className).lookup("path"));
            word classFile(parser.subDict(className).lookup("name"));

            Info<< "Showing documentation for type " << className << nl << endl;

            Info<< "Source file: " << classFolder.c_str() << classFile << nl
                << endl;

            system(docBrowser);
        }
        else
        {
            FatalErrorIn
            (
                "void Foam::helpType::displayDoc(const word, const string)"
            )
                << "No help for type " << className << " found."
                << "  Valid options include:"
                << SortableList<word>(parser.toc())
                << exit(FatalError);
        }
    }
    else
    {
        Info<< "No Doxygen sources found under search paths: "
            << docDirs << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helpType::helpType()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::helpType::~helpType()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::helpType::init()
{
    argList::addOption
    (
        "browse",
        "word",
        "display documentation for boundary condition in browser"
    );
}


// ************************************************************************* //

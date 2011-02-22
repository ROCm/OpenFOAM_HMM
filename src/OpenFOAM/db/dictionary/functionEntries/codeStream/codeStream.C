/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "codeStream.H"
#include "addToMemberFunctionSelectionTable.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "IOstreams.H"
#include "SHA1Digest.H"
#include "OSHA1stream.H"
#include "codeStreamTools.H"
#include "stringOps.H"
#include "dlLibraryTable.H"
#include "OSspecific.H"
#include "Time.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(codeStream, 0);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        codeStream,
        execute,
        primitiveEntryIstream
    );

}
}


const Foam::word Foam::functionEntries::codeStream::codeTemplateC
    = "codeStreamTemplate.C";


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::codeStream::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    if (isAdministrator())
    {
        FatalIOErrorIn
        (
            "functionEntries::codeStream::execute(..)",
            parentDict
        )   << "This code should not be executed by someone with administrator"
            << " rights due to security reasons." << endl
            << "(it writes a shared library which then gets loaded "
            << "using dlopen)"
            << exit(FatalIOError);
    }


    // Read three sections of code.
    // Remove any leading whitespace - necessary for compilation options,
    // convenience for includes and body.
    dictionary codeDict(is);

    // "codeInclude" is optional
    string codeInclude;
    if (codeDict.found("codeInclude"))
    {
        codeInclude = stringOps::trimLeft(codeDict["codeInclude"]);
    }

    // "codeOptions" is optional
    string codeOptions;
    if (codeDict.found("codeOptions"))
    {
        codeOptions = stringOps::trimLeft(codeDict["codeOptions"]);
    }

    // "code" is mandatory
    string code = stringOps::trimLeft(codeDict["code"]);

    // Create name from the contents
    SHA1Digest sha;
    {
        OSHA1stream os;
        os  << codeInclude << codeOptions << code;
        sha = os.digest();
    }

    // write code into _SHA1 subdirectory
    fileName dir;
    if (isA<IOdictionary>(parentDict))
    {
        const IOdictionary& d = static_cast<const IOdictionary&>(parentDict);
        dir = d.db().time().constantPath()/"codeStream"/"_" + sha.str();
    }
    else
    {
        dir = fileName("codeStream")/"_" + sha.str();
    }

    fileName name = "codeStream_" + sha.str();
    fileName libPath
    (
        Foam::getEnv("FOAM_USER_LIBBIN")/"lib" + name + ".so"
    );

    void* lib = dlLibraryTable::findLibrary(libPath);

    if (!lib)
    {
        if (Pstream::master())
        {
            if (!codeStreamTools::upToDate(dir, sha))
            {
                Info<< "Creating new library in " << libPath << endl;

                const fileName fileCsrc
                (
                    codeStreamTools::findTemplate
                    (
                        codeTemplateC
                    )
                );

                // not found!
                if (fileCsrc.empty())
                {
                    FatalIOErrorIn
                    (
                        "functionEntries::codeStream::execute(..)",
                        parentDict
                    )   << "Could not find the code template: "
                        << codeTemplateC << nl
                        << codeStreamTools::searchedLocations()
                        << exit(FatalIOError);
                }


                List<codeStreamTools::fileAndVars> copyFiles(1);
                copyFiles[0].file() = fileCsrc;
                copyFiles[0].set("codeInclude", codeInclude);
                copyFiles[0].set("code", code);

                List<codeStreamTools::fileAndContent> filesContents(2);

                // Write Make/files
                filesContents[0].first() = "Make/files";
                filesContents[0].second() =
                    codeTemplateC + "\n\n"
                    "LIB = $(FOAM_USER_LIBBIN)/lib" + name;

                // Write Make/options
                filesContents[1].first() = "Make/options";
                filesContents[1].second() =
                    "EXE_INC = -g \\\n" + codeOptions + "\n\nLIB_LIBS =";

                codeStreamTools writer(name, copyFiles, filesContents);
                if (!writer.copyFilesContents(dir))
                {
                    FatalIOErrorIn
                    (
                        "functionEntries::codeStream::execute(..)",
                        parentDict
                    )   << "Failed writing " <<nl
                        << copyFiles << endl
                        << filesContents
                        << exit(FatalIOError);
                }
            }

            Foam::string wmakeCmd("wmake libso " + dir);
            Info<< "Invoking " << wmakeCmd << endl;
            if (Foam::system(wmakeCmd))
            {
                FatalIOErrorIn
                (
                    "functionEntries::codeStream::execute(..)",
                    parentDict
                )   << "Failed " << wmakeCmd
                    << exit(FatalIOError);
            }
        }

        if (!dlLibraryTable::open(libPath))
        {
            FatalIOErrorIn
            (
                "functionEntries::codeStream::execute(..)",
                parentDict
            )   << "Failed loading library " << libPath
                << exit(FatalIOError);
        }

        lib = dlLibraryTable::findLibrary(libPath);
    }
    else
    {
        Info<< "Reusing library in " << libPath << endl;
    }


    // Find the library handle.
    void (*function)(const dictionary&, Ostream&);
    function = reinterpret_cast<void(*)(const dictionary&, Ostream&)>
    (
        dlSym(lib, name)
    );

    if (!function)
    {
        FatalIOErrorIn
        (
            "functionEntries::codeStream::execute(..)",
            parentDict
        )   << "Failed looking up symbol " << name
            << " in library " << lib << exit(FatalIOError);
    }

    OStringStream os;
    (*function)(parentDict, os);
    IStringStream resultStream(os.str());
    entry.read(parentDict, resultStream);

    return true;
}


// ************************************************************************* //

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

    // get code dictionary
    // must reference parent for stringOps::expand to work nicely
    dictionary codeDict("#codeStream", parentDict, is);


    // Read three sections of code.
    // Remove any leading whitespace - necessary for compilation options,
    // convenience for includes and body.

    // "codeInclude" is optional
    string codeInclude;
    if (codeDict.found("codeInclude"))
    {
        codeInclude = stringOps::trim(codeDict["codeInclude"]);
        stringOps::inplaceExpand(codeInclude, codeDict);
    }

    // "codeOptions" is optional
    string codeOptions;
    if (codeDict.found("codeOptions"))
    {
        codeOptions = stringOps::trim(codeDict["codeOptions"]);
        stringOps::inplaceExpand(codeOptions, codeDict);
    }

    // "code" is mandatory
    string code = stringOps::trim(codeDict["code"]);
    stringOps::inplaceExpand(code, codeDict);


    // Create SHA1 digest from the contents
    SHA1Digest sha;
    {
        OSHA1stream os;
        os  << codeInclude << codeOptions << code;
        sha = os.digest();
    }


    // codeName = prefix + sha1
    const fileName codeName = "codeStream_" + sha.str();

    // write code into _SHA1 subdir
    const fileName codePath = codeStreamTools::codePath("_" + sha.str());

    // write library into platforms/$WM_OPTIONS/lib subdir
    const fileName libPath = codeStreamTools::libPath(codeName);


    void* lib = dlLibraryTable::findLibrary(libPath);

    // try to load if not already loaded
    if (!lib && dlLibraryTable::open(libPath, false))
    {
        lib = dlLibraryTable::findLibrary(libPath);
    }

    // did not load - need to compile it
    if (!lib)
    {
        if (Pstream::master())
        {
            if (!codeStreamTools::upToDate(codePath, sha))
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
                  + codeStreamTools::libTarget(codeName);

                // Write Make/options
                filesContents[1].first() = "Make/options";
                filesContents[1].second() =
                    "EXE_INC = -g \\\n"
                  + codeOptions
                  + "\n\nLIB_LIBS =";

                codeStreamTools writer(codeName, copyFiles, filesContents);
                if (!writer.copyFilesContents(codePath))
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

            const Foam::string wmakeCmd("wmake libso " + codePath);
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

//        bool dummy = true;
//        reduce(dummy, orOp<bool>());

        if (!dlLibraryTable::open(libPath, false))
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
        dlSym(lib, codeName)
    );

    if (!function)
    {
        FatalIOErrorIn
        (
            "functionEntries::codeStream::execute(..)",
            parentDict
        )   << "Failed looking up symbol " << codeName
            << " in library " << lib << exit(FatalIOError);
    }

    OStringStream os(is.format());
    (*function)(parentDict, os);
    IStringStream resultStream(os.str());
    entry.read(parentDict, resultStream);

    return true;
}


// ************************************************************************* //

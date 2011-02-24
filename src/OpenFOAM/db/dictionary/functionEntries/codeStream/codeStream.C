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
#include "stringOps.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "dlLibraryTable.H"
#include "OSspecific.H"
#include "Time.H"
#include "PstreamReduceOps.H"


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
    dynamicCode::checkSecurity
    (
        "functionEntries::codeStream::execute(..)",
        parentDict
    );

    // get code dictionary
    // must reference parent for stringOps::expand to work nicely
    dictionary codeDict("#codeStream", parentDict, is);

    // get code, codeInclude, codeOptions
    dynamicCodeContext context(codeDict);

    // codeName: prefix_ + sha1
    // codeDir : _<sha1>
    dynamicCode dynCode
    (
        "codeStream_" + context.sha1().str(),
        "_" + context.sha1().str()
    );

    // Load library if not already loaded
    // Version information is encoded in the libPath (encoded with the SHA1)
    const fileName libPath = dynCode.libPath();

    void* lib = dlLibraryTable::findLibrary(libPath);

    if (!lib && dlLibraryTable::open(libPath, false))
    {
        lib = dlLibraryTable::findLibrary(libPath);
    }

    // did not load - need to compile it
    if (!lib)
    {
        if (Pstream::master())
        {
            if (!dynCode.upToDate(context))
            {
                Info<< "Creating new library in "
                    << dynCode.libPath() << endl;

                // filter C template
                dynCode.addFilterFile(codeTemplateC);

                // filter with this context
                dynCode.setFilterContext(context);

                // Write Make/files
                dynCode.addCreateFile
                (
                    "Make/files",
                    codeTemplateC + "\n\n"
                  + dynCode.libTarget()
                );

                // Write Make/options
                dynCode.addCreateFile
                (
                    "Make/options",
                    "EXE_INC = -g \\\n"
                  + context.options()
                  + "\n\nLIB_LIBS ="
                );

                if (!dynCode.copyFilesContents())
                {
                    FatalIOErrorIn
                    (
                        "functionEntries::codeStream::execute(..)",
                        parentDict
                    )   << "Failed writing " << nl
                        // << copyFiles << endl
                        // << filesContents
                        << exit(FatalIOError);
                }
            }

            if (!dynCode.wmakeLibso())
            {
                FatalIOErrorIn
                (
                    "functionEntries::codeStream::execute(..)",
                    parentDict
                )   << "Failed wmake " << libPath
                    << exit(FatalIOError);
            }
        }

        // all processes must wait for compile
        bool dummy = true;
        reduce(dummy, orOp<bool>());

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
        dlSym(lib, dynCode.codeName())
    );


    if (!function)
    {
        FatalIOErrorIn
        (
            "functionEntries::codeStream::execute(..)",
            parentDict
        )   << "Failed looking up symbol " << dynCode.codeName()
            << " in library " << lib << exit(FatalIOError);
    }

    OStringStream os(is.format());
    (*function)(parentDict, os);
    IStringStream resultStream(os.str());
    entry.read(parentDict, resultStream);

    return true;
}


// ************************************************************************* //

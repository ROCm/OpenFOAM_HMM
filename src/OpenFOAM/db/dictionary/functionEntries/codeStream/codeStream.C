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
    Info<< "Using #codeStream at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << endl;

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

    // codeName: codeStream + _<sha1>
    // codeDir : _<sha1>
    dynamicCode dynCode
    (
        "codeStream" + context.sha1().str(true),
        context.sha1().str(true)
    );

    // Load library if not already loaded
    // Version information is encoded in the libPath (encoded with the SHA1)
    const fileName libPath = dynCode.libPath();

    // see if library is loaded
    void* lib = dlLibraryTable::findLibrary(libPath);


    if (!lib)
    {
        Info<< "Using #codeStream with " << libPath << endl;
    }


    // nothing loaded
    // avoid compilation if possible by loading an existing library
    if (!lib && dlLibraryTable::open(libPath, false))
    {
        lib = dlLibraryTable::findLibrary(libPath);
    }


    // create library if required
    if (!lib)
    {
        bool create = Pstream::master();

        if (create)
        {
            if (!dynCode.upToDate(context))
            {
                // filter with this context
                dynCode.reset(context);

                // compile filtered C template
                dynCode.addCompileFile(codeTemplateC);

                // define Make/options
                dynCode.setMakeOptions
                (
                    "EXE_INC = -g \\\n"
                  + context.options()
                  + "\n\nLIB_LIBS ="
                );

                if (!dynCode.copyOrCreateFiles(true))
                {
                    FatalIOErrorIn
                    (
                        "functionEntries::codeStream::execute(..)",
                        parentDict
                    )   << "Failed writing files for" << nl
                        << dynCode.libRelPath() << nl
                        << exit(FatalIOError);
                }
            }

            if (!dynCode.wmakeLibso())
            {
                FatalIOErrorIn
                (
                    "functionEntries::codeStream::execute(..)",
                    parentDict
                )   << "Failed wmake " << dynCode.libRelPath() << nl
                    << exit(FatalIOError);
            }
        }

        // all processes must wait for compile to finish
        reduce(create, orOp<bool>());

        if (!dlLibraryTable::open(libPath, false))
        {
            FatalIOErrorIn
            (
                "functionEntries::codeStream::execute(..)",
                parentDict
            )   << "Failed loading library " << libPath << nl
                << exit(FatalIOError);
        }

        lib = dlLibraryTable::findLibrary(libPath);
    }


    // Find the function handle in the library
    streamingFunctionType function =
        reinterpret_cast<streamingFunctionType>
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

    // use function to write stream
    OStringStream os(is.format());
    (*function)(os, parentDict);

    // get the entry from this stream
    IStringStream resultStream(os.str());
    entry.read(parentDict, resultStream);

    return true;
}


// ************************************************************************* //

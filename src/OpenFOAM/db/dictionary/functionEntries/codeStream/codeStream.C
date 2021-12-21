/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "StringStream.H"
#include "Time.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(codeStream, 0);

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        codeStream,
        execute,
        dictionaryIstream,
        codeStream
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        codeStream,
        execute,
        primitiveEntryIstream,
        codeStream
    );
} // End namespace functionEntries
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dlLibraryTable& Foam::functionEntries::codeStream::libs
(
    const dictionary& dict
)
{
    return static_cast<const baseIOdictionary&>(dict.topDict()).time().libs();
}


bool Foam::functionEntries::codeStream::doingMasterOnlyReading
(
    const dictionary& dict
)
{
    // Fallback value
    bool masterOnly = regIOobject::masterOnlyReading;

    const auto* iodictPtr = isA<baseIOdictionary>(dict.topDict());

    if (iodictPtr)
    {
        masterOnly = iodictPtr->globalObject();

        DebugPout
            << "codeStream : baseIOdictionary:" << dict.name()
            << " master-only-reading:" << masterOnly << endl;
    }
    else
    {
        DebugPout
            << "codeStream : not a baseIOdictionary:" << dict.name()
            << " master-only-reading:" << masterOnly << endl;
    }

    return masterOnly;
}


Foam::functionEntries::codeStream::streamingFunctionType
Foam::functionEntries::codeStream::getFunction
(
    const dictionary& parentDict,
    const dictionary& codeDict
)
{
    // get code, codeInclude, codeOptions
    dynamicCodeContext context(codeDict);

    // codeName: codeStream + _<sha1>
    // codeDir : _<sha1>
    std::string sha1Str(context.sha1().str(true));
    dynamicCode dynCode("codeStream" + sha1Str, sha1Str);

    // Load library if not already loaded
    // Version information is encoded in the libPath (encoded with the SHA1)
    const fileName libPath = dynCode.libPath();

    // see if library is loaded
    void* lib = nullptr;

    const dictionary& topDict = parentDict.topDict();

    if (isA<baseIOdictionary>(topDict))
    {
        lib = libs(parentDict).findLibrary(libPath);
    }

    // nothing loaded
    // avoid compilation if possible by loading an existing library
    if (!lib)
    {
        DetailInfo
            << "Using #codeStream with " << libPath << endl;

        if (isA<baseIOdictionary>(topDict))
        {
            // Cached access to libs, with cleanup upon termination
            lib = libs(parentDict).open(libPath, false);
        }
        else
        {
            // Uncached opening of libPath. Do not complain if cannot be loaded
            lib = Foam::dlOpen(libPath, false);
        }
    }


    // create library if required
    if (!lib)
    {
        const bool create =
            Pstream::master()
         || (IOobject::fileModificationSkew <= 0);   // not NFS

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
                  + "\n\nLIB_LIBS = \\\n"
                    "    -lOpenFOAM \\\n"
                  + context.libs()
                );

                if (!dynCode.copyOrCreateFiles(true))
                {
                    FatalIOErrorInFunction(parentDict)
                        << "Failed writing files for" << nl
                        << dynCode.libRelPath() << nl
                        << exit(FatalIOError);
                }
            }

            if (!dynCode.wmakeLibso())
            {
                FatalIOErrorInFunction(parentDict)
                    << "Failed wmake " << dynCode.libRelPath() << nl
                    << exit(FatalIOError);
            }
        }

        //- Only block if we're not doing master-only reading. (flag set by
        //  regIOobject::read, baseIOdictionary constructor)
        if
        (
           !doingMasterOnlyReading(topDict)
         && IOobject::fileModificationSkew > 0
        )
        {
            //- Since the library has only been compiled on the master the
            //  other nodes need to pick this library up through NFS
            //  We do this by just polling a few times using the
            //  fileModificationSkew.

            off_t mySize = Foam::fileSize(libPath);
            off_t masterSize = mySize;
            Pstream::scatter(masterSize);

            for
            (
                label iter = 0;
                iter < IOobject::maxFileModificationPolls;
                ++iter
            )
            {
                DebugPout
                    << "on processor " << Pstream::myProcNo()
                    << "masterSize:" << masterSize
                    << " and localSize:" << mySize
                    << endl;

                if (mySize == masterSize)
                {
                    break;
                }
                else if (mySize > masterSize)
                {
                    FatalIOErrorInFunction(context.dict())
                        << "Excessive size when reading (NFS mounted) library "
                        << nl << libPath << nl
                        << "on processor " << Pstream::myProcNo()
                        << " detected size " << mySize
                        << " whereas master size is " << masterSize
                        << " bytes." << nl
                        << "If your case is NFS mounted increase"
                        << " fileModificationSkew or maxFileModificationPolls;"
                        << nl << "If your case is not NFS mounted"
                        << " (so distributed) set fileModificationSkew"
                        << " to 0"
                        << exit(FatalIOError);
                }
                else
                {
                    DebugPout
                        << "Local file " << libPath
                        << " not of same size (" << mySize
                        << ") as master ("
                        << masterSize << "). Waiting for "
                        << IOobject::fileModificationSkew
                        << " seconds." << endl;

                    Foam::sleep(IOobject::fileModificationSkew);

                    // Recheck local size
                    mySize = Foam::fileSize(libPath);
                }
            }


            // Finished doing iterations. Do final check
            if (mySize != masterSize)
            {
                FatalIOErrorInFunction(context.dict())
                    << "Cannot read (NFS mounted) library " << nl
                    << libPath << nl
                    << "on processor " << Pstream::myProcNo()
                    << " detected size " << mySize
                    << " whereas master size is " << masterSize
                    << " bytes." << nl
                    << "If your case is NFS mounted increase"
                    << " fileModificationSkew or maxFileModificationPolls;"
                    << nl << "If your case is not NFS mounted"
                    << " (so distributed) set fileModificationSkew"
                    << " to 0"
                    << exit(FatalIOError);
            }

            DebugPout
                << "on processor " << Pstream::myProcNo()
                << " after waiting: have masterSize:" << masterSize
                << " and localSize:" << mySize << endl;
        }

        if (isA<baseIOdictionary>(topDict))
        {
            // Cached access to libs, with cleanup upon termination
            DebugPout
                << "Opening cached dictionary:" << libPath << endl;

            lib = libs(parentDict).open(libPath, false);

            if (!lib)
            {
                FatalIOErrorInFunction(parentDict)
                    << "Failed loading library " << libPath << nl
                    << "Did you add all libraries to the 'libs' entry"
                    << " in system/controlDict?"
                    << exit(FatalIOError);
            }
        }
        else
        {
            // Uncached opening of libPath
            DebugPout
                << "Opening uncached dictionary:" << libPath << endl;

            lib = Foam::dlOpen(libPath, true);
        }
    }

    bool haveLib = lib;
    if (!doingMasterOnlyReading(topDict))
    {
        reduce(haveLib, andOp<bool>());
    }

    if (!haveLib)
    {
        FatalIOErrorInFunction(parentDict)
            << "Failed loading library " << libPath
            << " on some processors."
            << exit(FatalIOError);
    }


    // Find the function handle in the library
    streamingFunctionType function =
        reinterpret_cast<streamingFunctionType>
        (
            Foam::dlSym(lib, dynCode.codeName())
        );


    if (!function)
    {
        FatalIOErrorInFunction(parentDict)
            << "Failed looking up symbol " << dynCode.codeName()
            << " in library " << lib << exit(FatalIOError);
    }

    return function;
}


Foam::string Foam::functionEntries::codeStream::evaluate
(
    const dictionary& parentDict,
    Istream& is
)
{
    DetailInfo
        << "Using #codeStream at line " << is.lineNumber()
        << " in file " <<  parentDict.relativeName() << endl;

    dynamicCode::checkSecurity
    (
        "functionEntries::codeStream::evaluate(..)",
        parentDict
    );

    // Get code dictionary
    dictionary codeDict("#codeStream", parentDict, is);

    // Use function to write stream
    OStringStream os(is.format());

    streamingFunctionType function = getFunction(parentDict, codeDict);
    (*function)(os, parentDict);

    // Return evaluated content as string
    return os.str();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::codeStream::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    IStringStream result(evaluate(parentDict, is));
    entry.read(parentDict, result);

    return true;
}


bool Foam::functionEntries::codeStream::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    IStringStream result(evaluate(parentDict, is));
    parentDict.read(result);

    return true;
}


// ************************************************************************* //

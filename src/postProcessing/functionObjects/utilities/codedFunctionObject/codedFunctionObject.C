/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "codedFunctionObject.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "SHA1Digest.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "stringOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(codedFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        codedFunctionObject,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void* Foam::codedFunctionObject::loadLibrary
(
    const fileName& libPath,
    const string& globalFuncName,
    const dictionary& contextDict
) const
{
    void* lib = 0;

    // avoid compilation by loading an existing library
    if (!libPath.empty())
    {
        dlLibraryTable& libs = const_cast<Time&>(time_).libs();

        if (libs.open(libPath, false))
        {
            lib = libs.findLibrary(libPath);

            // verify the loaded version and unload if needed
            if (lib)
            {
                // provision for manual execution of code after loading
                if (dlSymFound(lib, globalFuncName))
                {
                    loaderFunctionType function =
                        reinterpret_cast<loaderFunctionType>
                        (
                            dlSym(lib, globalFuncName)
                        );

                    if (function)
                    {
                        (*function)(true);    // force load
                    }
                    else
                    {
                        FatalIOErrorIn
                        (
                            "codedFunctionObject::updateLibrary()",
                            contextDict
                        )   << "Failed looking up symbol " << globalFuncName
                            << nl << "from " << libPath << exit(FatalIOError);
                    }
                }
                else
                {
                    FatalIOErrorIn
                    (
                        "codedFunctionObject::loadLibrary()",
                        contextDict
                    )   << "Failed looking up symbol " << globalFuncName << nl
                        << "from " << libPath << exit(FatalIOError);

                    lib = 0;
                    if (!libs.close(libPath, false))
                    {
                        FatalIOErrorIn
                        (
                            "codedFunctionObject::loadLibrary()",
                            contextDict
                        )   << "Failed unloading library "
                            << libPath
                            << exit(FatalIOError);
                    }
                }
            }
        }
    }

    return lib;
}


void Foam::codedFunctionObject::unloadLibrary
(
    const fileName& libPath,
    const string& globalFuncName,
    const dictionary& contextDict
) const
{
    void* lib = 0;

    if (libPath.empty())
    {
        return;
    }

    dlLibraryTable& libs = const_cast<Time&>(time_).libs();

    lib = libs.findLibrary(libPath);

    if (!lib)
    {
        return;
    }

    // provision for manual execution of code before unloading
    if (dlSymFound(lib, globalFuncName))
    {
        loaderFunctionType function =
            reinterpret_cast<loaderFunctionType>
            (
                dlSym(lib, globalFuncName)
            );

        if (function)
        {
            (*function)(false);    // force unload
        }
        else
        {
            FatalIOErrorIn
            (
                "codedFunctionObject::unloadLibrary()",
                contextDict
            )   << "Failed looking up symbol " << globalFuncName << nl
                << "from " << libPath << exit(FatalIOError);
        }
    }

    if (!libs.close(libPath, false))
    {
        FatalIOErrorIn
        (
            "codedFunctionObject::"
            "updateLibrary()",
            contextDict
        )   << "Failed unloading library " << libPath
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::codedFunctionObject::createLibrary
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    bool create = Pstream::master();

    if (create)
    {
        // Write files for new library
        if (!dynCode.upToDate(context))
        {
            // filter with this context
            dynCode.reset(context);

            // Set additional rewrite rules
            dynCode.setFilterVariable("typeName", redirectType_);
            dynCode.setFilterVariable("codeRead", codeRead_);
            dynCode.setFilterVariable("codeExecute", codeExecute_);
            dynCode.setFilterVariable("codeEnd", codeEnd_);
            //dynCode.setFilterVariable("codeWrite", codeWrite_);

            // compile filtered C template
            dynCode.addCompileFile("functionObjectTemplate.C");
            dynCode.addCompileFile("FilterFunctionObjectTemplate.C");

            // copy filtered H template
            dynCode.addCopyFile("FilterFunctionObjectTemplate.H");
            dynCode.addCopyFile("functionObjectTemplate.H");
            dynCode.addCopyFile("IOfunctionObjectTemplate.H");

            // debugging: make BC verbose
            //         dynCode.setFilterVariable("verbose", "true");
            //         Info<<"compile " << redirectType_ << " sha1: "
            //             << context.sha1() << endl;

            // define Make/options
            dynCode.setMakeOptions
            (
                "EXE_INC = -g \\\n"
                "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
              + context.options()
              + "\n\nLIB_LIBS = \\\n"
              + "    -lOpenFOAM \\\n"
              + "    -lfiniteVolume \\\n"
              + context.libs()
            );

            if (!dynCode.copyOrCreateFiles(true))
            {
                FatalIOErrorIn
                (
                    "codedFunctionObject::createLibrary(..)",
                    context.dict()
                )   << "Failed writing files for" << nl
                    << dynCode.libRelPath() << nl
                    << exit(FatalIOError);
            }
        }

        if (!dynCode.wmakeLibso())
        {
            FatalIOErrorIn
            (
                "codedFunctionObject::createLibrary(..)",
                context.dict()
            )   << "Failed wmake " << dynCode.libRelPath() << nl
                << exit(FatalIOError);
        }
    }


    // all processes must wait for compile to finish
    reduce(create, orOp<bool>());
}


void Foam::codedFunctionObject::updateLibrary() const
{
    dynamicCode::checkSecurity
    (
        "codedFunctionObject::updateLibrary()",
        dict_
    );

    dynamicCodeContext context(dict_);

    // codeName: redirectType + _<sha1>
    // codeDir : redirectType
    dynamicCode dynCode
    (
        redirectType_ + context.sha1().str(true),
        redirectType_
    );
    const fileName libPath = dynCode.libPath();


    // the correct library was already loaded => we are done
    if (const_cast<Time&>(time_).libs().findLibrary(libPath))
    {
        return;
    }

    Info<< "Using dynamicCode for functionObject " << name()
        << " at line " << dict_.startLineNumber()
        << " in " << dict_.name() << endl;


    // remove instantiation of fvPatchField provided by library
    redirectFunctionObjectPtr_.clear();

    // may need to unload old library
    unloadLibrary
    (
        oldLibPath_,
        dynamicCode::libraryBaseName(oldLibPath_),
        context.dict()
    );

    // try loading an existing library (avoid compilation when possible)
    if (!loadLibrary(libPath, dynCode.codeName(), context.dict()))
    {
        createLibrary(dynCode, context);

        loadLibrary(libPath, dynCode.codeName(), context.dict());
    }

    // retain for future reference
    oldLibPath_ = libPath;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedFunctionObject::codedFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    functionObject(name),
    time_(time),
    dict_(dict)
{
    read(dict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::codedFunctionObject::~codedFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::functionObject&
Foam::codedFunctionObject::redirectFunctionObject() const
{
    if (!redirectFunctionObjectPtr_.valid())
    {
        dictionary constructDict(dict_);
        constructDict.set("type", redirectType_);

        redirectFunctionObjectPtr_ = functionObject::New
        (
            redirectType_,
            time_,
            constructDict
        );
    }
    return redirectFunctionObjectPtr_();
}


bool Foam::codedFunctionObject::start()
{
    updateLibrary();
    return redirectFunctionObject().start();
}


bool Foam::codedFunctionObject::execute(const bool forceWrite)
{
    updateLibrary();
    return redirectFunctionObject().execute(forceWrite);
}


bool Foam::codedFunctionObject::end()
{
    updateLibrary();
    return redirectFunctionObject().end();
}


bool Foam::codedFunctionObject::read(const dictionary& dict)
{
    dict.lookup("redirectType") >> redirectType_;

    const entry* readPtr = dict.lookupEntryPtr
    (
        "codeRead",
        false,
        false
    );
    if (readPtr)
    {
        codeRead_ = stringOps::trim(readPtr->stream());
        stringOps::inplaceExpand(codeRead_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeRead_,
            readPtr->startLineNumber(),
            dict.name()
        );
    }

    const entry* execPtr = dict.lookupEntryPtr
    (
        "codeExecute",
        false,
        false
    );
    if (execPtr)
    {
        codeExecute_ = stringOps::trim(execPtr->stream());
        stringOps::inplaceExpand(codeExecute_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeExecute_,
            execPtr->startLineNumber(),
            dict.name()
        );
    }

    const entry* endPtr = dict.lookupEntryPtr
    (
        "codeEnd",
        false,
        false
    );
    if (execPtr)
    {
        codeEnd_ = stringOps::trim(endPtr->stream());
        stringOps::inplaceExpand(codeEnd_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeEnd_,
            endPtr->startLineNumber(),
            dict.name()
        );
    }

    updateLibrary();
    return redirectFunctionObject().read(dict);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "codedFixedValueFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "dlLibraryTable.H"
#include "IFstream.H"
#include "OFstream.H"
#include "SHA1Digest.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::codedFixedValueFvPatchScalarField::codeTemplateC
    = "fixedValueFvPatchScalarFieldTemplate.C";

const Foam::word Foam::codedFixedValueFvPatchScalarField::codeTemplateH
    = "fixedValueFvPatchScalarFieldTemplate.H";

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


Foam::string Foam::codedFixedValueFvPatchScalarField::libraryGlobalName
(
    const fileName& libPath
)
{
    // global function name (SHA1-encoded)
    // that can be used for version control and/or explicit loading/unloading
    string globalFuncName = libPath.name();

    // remove ".so" extension
    string::size_type dot = globalFuncName.find('.');
    if (dot != string::npos)
    {
        globalFuncName.resize(dot);
    }

    // remove leading 'lib' from name
    globalFuncName.erase(0, 3);

    return globalFuncName;
}




void* Foam::codedFixedValueFvPatchScalarField::loadLibrary
(
    const fileName& libPath,
    const string& globalFuncName,
    const dictionary& contextDict
)
{
    void* lib = 0;

    // global function name (SHA1-encoded)
    // that can be used for version control and/or explicit loading/unloading

    // avoid compilation by loading an existing library
    if (dlLibraryTable::open(libPath, false))
    {
        lib = dlLibraryTable::findLibrary(libPath);

        // verify the loaded version and unload if needed
        if (lib)
        {
            if (dlSymFound(lib, globalFuncName))
            {
                // Find the function handle in the library
                loaderFunctionType function =
                    reinterpret_cast<loaderFunctionType>
                    (
                        dlSym(lib, globalFuncName)
                    );

                if (function)
                {
                    // force load
                    (*function)(true);
                }
                else
                {
                    FatalIOErrorIn
                    (
                        "codedFixedValueFvPatchScalarField::updateLibrary()",
                        contextDict
                    )   << "Failed looking up symbol " << globalFuncName << nl
                        << "from " << libPath << exit(FatalIOError);
                }
            }
            else
            {
                FatalIOErrorIn
                (
                    "codedFixedValueFvPatchScalarField::loadLibrary()",
                    contextDict
                )   << "Failed looking up symbol " << globalFuncName << nl
                    << "from " << libPath << exit(FatalIOError);

                lib = 0;
                if (!dlLibraryTable::close(libPath, false))
                {
                    FatalIOErrorIn
                    (
                        "codedFixedValueFvPatchScalarField::loadLibrary()",
                        contextDict
                    )   << "Failed unloading library "
                        << libPath
                        << exit(FatalIOError);
                }
            }
        }
    }

    return lib;
}


void Foam::codedFixedValueFvPatchScalarField::unloadLibrary
(
    const fileName& libPath,
    const string& globalFuncName,
    const dictionary& contextDict
)
{
    void* lib = dlLibraryTable::findLibrary(libPath);

    if (!lib)
    {
        return;
    }

    // provision for manual execution of code before unloading
    if (dlSymFound(lib, globalFuncName))
    {
        // Find the function handle in the library
        loaderFunctionType function =
            reinterpret_cast<loaderFunctionType>
            (
                dlSym(lib, globalFuncName)
            );

        if (function)
        {
            // force unload
            (*function)(false);
        }
        else
        {
            FatalIOErrorIn
            (
                "codedFixedValueFvPatchScalarField::unloadLibrary()",
                contextDict
            )   << "Failed looking up symbol " << globalFuncName << nl
                << "from " << libPath
                << exit(FatalIOError);
        }
    }

    if (!dlLibraryTable::close(libPath, false))
    {
        FatalIOErrorIn
        (
            "codedFixedValueFvPatchScalarField::"
            "updateLibrary()",
            contextDict
        )   << "Failed unloading library " << libPath
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::IOdictionary& Foam::codedFixedValueFvPatchScalarField::dict() const
{
    if (db().foundObject<IOdictionary>("codeDict"))
    {
        return db().lookupObject<IOdictionary>("codeDict");
    }
    else
    {
        return db().store
        (
            new IOdictionary
            (
                IOobject
                (
                    "codeDict",
                    db().time().system(),
                    db(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


void Foam::codedFixedValueFvPatchScalarField::createLibrary
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Write files for new library
    if (Pstream::master() && !dynCode.upToDate(context))
    {
        // filter with this context
        dynCode.reset(context);

        // compile filtered C template
        dynCode.addCompileFile(codeTemplateC);

        // copy filtered H template
        dynCode.addCopyFile(codeTemplateH);

        // take no chances - typeName must be identical to redirectType_
        dynCode.setFilterVariable("typeName", redirectType_);

        // debugging: make BC verbose
//         dynCode.setFilterVariable("verbose", "true");
//         Info<<"compile " << redirectType_ << " sha1: "
//             << context.sha1() << endl;

        // define Make/options
        dynCode.setMakeOptions
        (
            "EXE_INC = -g \\\n"
            "-I$(LIB_SRC)/finiteVolume/lnInclude\\\n"
          + context.options()
          + "\n\nLIB_LIBS = "
        );

        if (!dynCode.copyOrCreateFiles(true))
        {
            FatalIOErrorIn
            (
                "codedFixedValueFvPatchScalarField::writeLibrary(..)",
                context.dict()
            )   << "Failed writing files for" << nl
                << dynCode.libPath() << nl
                << exit(FatalIOError);
        }
    }
}


void Foam::codedFixedValueFvPatchScalarField::updateLibrary() const
{
    dynamicCode::checkSecurity
    (
        "codedFixedValueFvPatchScalarField::updateLibrary()",
        dict_
    );

    // use system/codeDict or in-line
    const dictionary& codeDict =
    (
        dict_.found("code")
      ? dict_
      : this->dict().subDict(redirectType_)
    );

    dynamicCodeContext context(codeDict);


    // codeName: redirectType + _<sha1>
    // codeDir : redirectType
    dynamicCode dynCode
    (
        redirectType_ + context.sha1().str(true),
        redirectType_
    );
    const fileName libPath = dynCode.libPath();


    // see if library is loaded
    void* lib = dlLibraryTable::findLibrary(libPath);

    // library not loaded, and also need to unload old version
    if (!lib && !oldLibPath_.empty())
    {
        // unload code
        // Remove instantiation of fvPatchField provided by library
        redirectPatchFieldPtr_.clear();

        unloadLibrary
        (
            oldLibPath_,
            libraryGlobalName(oldLibPath_),
            context.dict()
        );
    }


    // retain for future reference
    oldLibPath_ = libPath;

    bool waiting = false;

    if (lib)
    {
        return;
    }

    // Remove instantiation of fvPatchField provided by library
    redirectPatchFieldPtr_.clear();

    // avoid compilation by loading an existing library
    lib = loadLibrary
    (
        libPath,
        dynCode.codeName(),
        context.dict()
    );


    // really do need to create library
    if (!lib)
    {
        if (Pstream::master())
        {
            createLibrary(dynCode, context);

            if (!dynCode.wmakeLibso())
            {
                FatalIOErrorIn
                (
                    "codedFixedValueFvPatchScalarField::updateLibrary()",
                    context.dict()
                )   << "Failed wmake " << libPath
                    << exit(FatalIOError);
            }
        }

        // all processes must wait for compile
        waiting = true;
        reduce(waiting, orOp<bool>());

        // load newly compiled libray
        lib = loadLibrary
        (
            libPath,
            dynCode.codeName(),
            context.dict()
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedFixedValueFvPatchScalarField::
codedFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    oldLibPath_(),
    redirectPatchFieldPtr_()
{}


Foam::codedFixedValueFvPatchScalarField::
codedFixedValueFvPatchScalarField
(
    const codedFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
    oldLibPath_(),
    redirectPatchFieldPtr_()
{}


Foam::codedFixedValueFvPatchScalarField::
codedFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    dict_(dict),
    redirectType_(dict.lookup("redirectType")),
    oldLibPath_(),
    redirectPatchFieldPtr_()
{
    updateLibrary();
}


Foam::codedFixedValueFvPatchScalarField::
codedFixedValueFvPatchScalarField
(
    const codedFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
    oldLibPath_(),
    redirectPatchFieldPtr_()
{}


Foam::codedFixedValueFvPatchScalarField::
codedFixedValueFvPatchScalarField
(
    const codedFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
    oldLibPath_(),
    redirectPatchFieldPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvPatchScalarField&
Foam::codedFixedValueFvPatchScalarField::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        // Construct a patch
        // Make sure to construct the patchfield with uptodate value.

        OStringStream os;
        os.writeKeyword("type") << redirectType_ << token::END_STATEMENT
            << nl;
        static_cast<const scalarField&>(*this).writeEntry("value", os);
        IStringStream is(os.str());
        dictionary dict(is);
//        Info<< "constructing patchField from :" << dict << endl;

//        if (fvPatchScalarField::dictionaryConstructorTablePtr_)
//        {
//            fvPatchScalarField::dictionaryConstructorPtr funcPtr =
//            (
//                fvPatchScalarField::dictionaryConstructorTablePtr_->
//                find(redirectType_)()
//            );
//
//            Info<< redirectType_ << " FunctionPtr => "
//                << long(funcPtr) << endl;
//        }

        redirectPatchFieldPtr_.set
        (
            fvPatchScalarField::New
            (
                patch(),
                dimensionedInternalField(),
                dict
            ).ptr()
        );
    }
    return redirectPatchFieldPtr_();
}


void Foam::codedFixedValueFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Make sure library containing user-defined fvPatchField is uptodate
    updateLibrary();

    const fvPatchScalarField& fvp = redirectPatchField();

    const_cast<fvPatchScalarField&>(fvp).updateCoeffs();

    // Copy through value
    operator==(fvp);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::codedFixedValueFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Make sure library containing user-defined fvPatchField is uptodate
    updateLibrary();

    const fvPatchScalarField& fvp = redirectPatchField();

    const_cast<fvPatchScalarField&>(fvp).evaluate(commsType);

    fixedValueFvPatchField<scalar>::evaluate(commsType);
}


void Foam::codedFixedValueFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
    os.writeKeyword("redirectType") << redirectType_
        << token::END_STATEMENT << nl;

    if (dict_.found("code"))
    {
        os.writeKeyword("code")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["code"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        codedFixedValueFvPatchScalarField
    );
}

// ************************************************************************* //

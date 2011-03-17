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
#include <cstring>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::codedFixedValueFvPatchScalarField::codeTemplateC
    = "fixedValueFvPatchFieldTemplate.C";

const Foam::word Foam::codedFixedValueFvPatchScalarField::codeTemplateH
    = "fixedValueFvPatchFieldTemplate.H";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void* Foam::codedFixedValueFvPatchScalarField::loadLibrary
(
    const fileName& libPath,
    const string& globalFuncName,
    const dictionary& contextDict
)
{
    void* lib = 0;

    // avoid compilation by loading an existing library
    if (!libPath.empty() && dlLibraryTable::open(libPath, false))
    {
        lib = dlLibraryTable::findLibrary(libPath);

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
    void* lib = 0;

    if (!libPath.empty())
    {
        lib = dlLibraryTable::findLibrary(libPath);
    }

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
                "codedFixedValueFvPatchScalarField::unloadLibrary()",
                contextDict
            )   << "Failed looking up symbol " << globalFuncName << nl
                << "from " << libPath << exit(FatalIOError);
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



template<class Type>
void Foam::codedFixedValueFvPatchScalarField::setFieldTemplates
(
    dynamicCode& dynCode
)
{
    word fieldType(pTraits<Type>::typeName);

    // template type for fvPatchField
    dynCode.setFilterVariable("TemplateType", fieldType);

    // Name for fvPatchField - eg, ScalarField, VectorField, ...
    fieldType[0] = toupper(fieldType[0]);
    dynCode.setFilterVariable("FieldType", fieldType + "Field");
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
    bool create = Pstream::master();

    if (create)
    {
        // Write files for new library
        if (!dynCode.upToDate(context))
        {
            // filter with this context
            dynCode.reset(context);

            // take no chances - typeName must be identical to redirectType_
            dynCode.setFilterVariable("typeName", redirectType_);

            // set TemplateType and FieldType filter variables
            // (for fvPatchField)
            setFieldTemplates<scalar>(dynCode);

            // compile filtered C template
            dynCode.addCompileFile(codeTemplateC);

            // copy filtered H template
            dynCode.addCopyFile(codeTemplateH);


            // debugging: make BC verbose
            //  dynCode.setFilterVariable("verbose", "true");
            //  Info<<"compile " << redirectType_ << " sha1: "
            //      << context.sha1() << endl;

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
                    "codedFixedValueFvPatchScalarField::createLibrary(..)",
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
                "codedFixedValueFvPatchScalarField::createLibrary(..)",
                context.dict()
            )   << "Failed wmake " << dynCode.libRelPath() << nl
                << exit(FatalIOError);
        }
    }


    // all processes must wait for compile to finish
    reduce(create, orOp<bool>());
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


    // the correct library was already loaded => we are done
    if (dlLibraryTable::findLibrary(libPath))
    {
        return;
    }

    // remove instantiation of fvPatchField provided by library
    redirectPatchFieldPtr_.clear();

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
    oldLibPath_(ptf.oldLibPath_),
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
    oldLibPath_(ptf.oldLibPath_),
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
    oldLibPath_(ptf.oldLibPath_),
    redirectPatchFieldPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvPatchScalarField&
Foam::codedFixedValueFvPatchScalarField::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        // Construct a patch
        // Make sure to construct the patchfield with up-to-date value

        OStringStream os;
        os.writeKeyword("type") << redirectType_ << token::END_STATEMENT
            << nl;
        static_cast<const scalarField&>(*this).writeEntry("value", os);
        IStringStream is(os.str());
        dictionary dict(is);

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

    // Make sure library containing user-defined fvPatchField is up-to-date
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
    // Make sure library containing user-defined fvPatchField is up-to-date
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

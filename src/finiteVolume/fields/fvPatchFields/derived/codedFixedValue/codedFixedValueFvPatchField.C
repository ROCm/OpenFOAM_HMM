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

#include "codedFixedValueFvPatchField.H"
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
#include "IOdictionary.H"

#include <dlfcn.h>
#include <link.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::word Foam::codedFixedValueFvPatchField<Type>::codeTemplateC
    = "fixedValueFvPatchFieldTemplate.C";

template<class Type>
const Foam::word Foam::codedFixedValueFvPatchField<Type>::codeTemplateH
    = "fixedValueFvPatchFieldTemplate.H";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
void* Foam::codedFixedValueFvPatchField<Type>::loadLibrary
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
        dlLibraryTable& libs = const_cast<Time&>(this->db().time()).libs();

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
                            "codedFixedValueFvPatchField<Type>::"
                            "updateLibrary()",
                            contextDict
                        )   << "Failed looking up symbol " << globalFuncName
                            << nl << "from " << libPath << exit(FatalIOError);
                    }
                }
                else
                {
                    FatalIOErrorIn
                    (
                        "codedFixedValueFvPatchField<Type>::loadLibrary()",
                        contextDict
                    )   << "Failed looking up symbol " << globalFuncName << nl
                        << "from " << libPath << exit(FatalIOError);

                    lib = 0;
                    if (!libs.close(libPath, false))
                    {
                        FatalIOErrorIn
                        (
                            "codedFixedValueFvPatchField<Type>::loadLibrary()",
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


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::unloadLibrary
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

    dlLibraryTable& libs = const_cast<Time&>(this->db().time()).libs();

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
                "codedFixedValueFvPatchField<Type>::unloadLibrary()",
                contextDict
            )   << "Failed looking up symbol " << globalFuncName << nl
                << "from " << libPath << exit(FatalIOError);
        }
    }

    if (!libs.close(libPath, false))
    {
        FatalIOErrorIn
        (
            "codedFixedValueFvPatchField<Type>::"
            "updateLibrary()",
            contextDict
        )   << "Failed unloading library " << libPath
            << exit(FatalIOError);
    }
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::setFieldTemplates
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

template<class Type>
const Foam::IOdictionary& Foam::codedFixedValueFvPatchField<Type>::dict() const
{
    const objectRegistry& obr = this->db();

    if (obr.foundObject<IOdictionary>("codeDict"))
    {
        return obr.lookupObject<IOdictionary>("codeDict");
    }
    else
    {
        return obr.store
        (
            new IOdictionary
            (
                IOobject
                (
                    "codeDict",
                    this->db().time().system(),
                    this->db(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::createLibrary
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
            setFieldTemplates(dynCode);

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
                    "codedFixedValueFvPatchField<Type>::createLibrary(..)",
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
                "codedFixedValueFvPatchField<Type>::createLibrary(..)",
                context.dict()
            )   << "Failed wmake " << dynCode.libRelPath() << nl
                << exit(FatalIOError);
        }
    }


    // all processes must wait for compile to finish
    reduce(create, orOp<bool>());
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::updateLibrary() const
{
    dynamicCode::checkSecurity
    (
        "codedFixedValueFvPatchField<Type>::updateLibrary()",
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
    if (const_cast<Time&>(this->db().time()).libs().findLibrary(libPath))
    {
        return;
    }

    Info<< "Using dynamicCode for patch " << this->patch().name()
        << " on field " << this->dimensionedInternalField().name() << nl
        << "at line " << codeDict.startLineNumber()
        << " in " << codeDict.name() << endl;


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

template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    oldLibPath_(),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const codedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
    oldLibPath_(ptf.oldLibPath_),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    dict_(dict),
    redirectType_(dict.lookup("redirectType")),
    oldLibPath_(),
    redirectPatchFieldPtr_()
{
    updateLibrary();
}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const codedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
    oldLibPath_(ptf.oldLibPath_),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const codedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
    oldLibPath_(ptf.oldLibPath_),
    redirectPatchFieldPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::fvPatchField<Type>&
Foam::codedFixedValueFvPatchField<Type>::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        // Construct a patch
        // Make sure to construct the patchfield with up-to-date value

        OStringStream os;
        os.writeKeyword("type") << redirectType_ << token::END_STATEMENT
            << nl;
        static_cast<const Field<Type>&>(*this).writeEntry("value", os);
        IStringStream is(os.str());
        dictionary dict(is);

        redirectPatchFieldPtr_.set
        (
            fvPatchField<Type>::New
            (
                this->patch(),
                this->dimensionedInternalField(),
                dict
            ).ptr()
        );
    }
    return redirectPatchFieldPtr_();
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary();

    const fvPatchField<Type>& fvp = redirectPatchField();

    const_cast<fvPatchField<Type>&>(fvp).updateCoeffs();

    // Copy through value
    this->operator==(fvp);

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary();

    const fvPatchField<Type>& fvp = redirectPatchField();

    const_cast<fvPatchField<Type>&>(fvp).evaluate(commsType);

    fixedValueFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fixedValueFvPatchField<Type>::write(os);
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


// ************************************************************************* //

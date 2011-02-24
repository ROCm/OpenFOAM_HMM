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
#include "codeProperties.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::codedFixedValueFvPatchScalarField::codeTemplateC
    = "fixedValueFvPatchScalarFieldTemplate.C";

const Foam::word Foam::codedFixedValueFvPatchScalarField::codeTemplateH
    = "fixedValueFvPatchScalarFieldTemplate.H";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::codeProperties&
Foam::codedFixedValueFvPatchScalarField::dict() const
{
    if (db().foundObject<codeProperties>(codeProperties::typeName))
    {
        return db().lookupObject<codeProperties>
        (
            codeProperties::typeName
        );
    }
    else
    {
        codeProperties* props = new codeProperties
        (
            IOobject
            (
                codeProperties::typeName,
                db().time().system(),
                db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

        return db().store(props);
    }
}


void Foam::codedFixedValueFvPatchScalarField::createLibrary
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
)
{
    // Write files for new library
    if (Pstream::master() && !dynCode.upToDate(context))
    {
        Info<< "Creating new library in " << dynCode.libPath() << endl;
        dynCode.clear();

        // filter with this context
        dynCode.setFilterContext(context);

        // filter C/H template
        dynCode.addFilterFile(codeTemplateC);
        dynCode.addFilterFile(codeTemplateH);

        // Make/files
        dynCode.addCreateFile
        (
            "Make/files",
            codeTemplateC + "\n\n"
          + dynCode.libTarget()
        );

        // Make/options
        dynCode.addCreateFile
        (
            "Make/options",
            "EXE_INC = -g \\\n"
            "-I$(LIB_SRC)/finiteVolume/lnInclude\\\n"
          + context.options()
          + "\n\nLIB_LIBS = "
        );

        if (!dynCode.copyFilesContents())
        {
            FatalIOErrorIn
            (
                "codedFixedValueFvPatchScalarField::writeLibrary(..)",
                context.dict()
            )   << "Failed writing " << nl
                // << copyFiles << nl
                // << filesContents
                << exit(FatalIOError);
        }
    }
}


void Foam::codedFixedValueFvPatchScalarField::updateLibrary()
{
    dynamicCode::checkSecurity
    (
        "codedFixedValueFvPatchScalarField::updateLibrary()",
        dict_
    );

    // use in-line or via codeProperties
    const bool useInlineDict = dict_.found("code");

    // determine code context (code, codeInclude, codeOptions)
    dynamicCodeContext context
    (
        useInlineDict
      ? dict_
      : this->dict().subDict(redirectType_)
    );


    // NOTE: probably don't need codeProperties anymore
    // since we use the sha1 directly
    if (!useInlineDict)
    {
        this->dict().setUnmodified();
    }


    // write code into redirectType_ subdir as well
    dynamicCode dynCode(redirectType_);

    // The version function name - based on the SHA1
    const string checkFuncName
    (
        dynCode.codeName() + "_" + context.sha1().str()
    );


    const fileName libPath = dynCode.libPath();
    void* lib = dlLibraryTable::findLibrary(libPath);

    // Load library if not already loaded
    bool reusing = false;
    if (!lib && dlLibraryTable::open(libPath, false))
    {
        lib = dlLibraryTable::findLibrary(libPath);
        reusing = true;
    }


    // library may have loaded, the version may not be correct
    bool waiting = false;
    if (lib)
    {
        // Unload library if needed
        if (!dlSym(lib, checkFuncName))
        {
            reusing = false;
            waiting = true;
            if (!dlLibraryTable::close(libPath, false))
            {
                FatalIOErrorIn
                (
                    "codedFixedValueFvPatchScalarField::updateLibrary(..)",
                    context.dict()
                )   << "Failed unloading library "
                    << libPath
                    << exit(FatalIOError);
            }
            lib = 0;
        }

        // unload from all processes
        reduce(waiting, orOp<bool>());
    }


    // create library
    if (!lib)
    {
        if (useInlineDict)
        {
            createLibrary(dynCode, context);
        }
        else
        {
            // Remove instantiation of fvPatchField provided by library
            redirectPatchFieldPtr_.clear();
            createLibrary(dynCode, context);
        }

        if (Pstream::master())
        {
            if (!dynCode.wmakeLibso())
            {
                FatalIOErrorIn
                (
                    "codedFixedValueFvPatchScalarField::updateLibrary()",
                    dict_
                )   << "Failed wmake " << libPath
                    << exit(FatalIOError);
            }
        }

        // all processes must wait for compile
        waiting = true;
        reduce(waiting, orOp<bool>());

        if (!dlLibraryTable::open(libPath, false))
        {
            FatalIOErrorIn
            (
                "codedFixedValueFvPatchScalarField::updateLibrary()",
                dict_
            )   << "Failed loading library " << libPath
                << exit(FatalIOError);
        }


        // paranoid - check that signature function is really there
        lib = dlLibraryTable::findLibrary(libPath);
        if (lib)
        {
            if (!dlSym(lib, checkFuncName))
            {
               FatalIOErrorIn
               (
                   "codedFixedValueFvPatchScalarField::updateLibrary(..)",
                   dict_
               )   << "Library loaded - but wrong version!"
                    << libPath
                    << exit(FatalIOError);
            }
            lib = 0;
        }
    }
    else if (reusing)
    {
        Info<< "Reusing library in " << libPath << nl
            << "    with " << context.sha1().str() << nl;
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
    redirectPatchFieldPtr_(NULL)
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
    redirectPatchFieldPtr_(NULL)
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
    redirectPatchFieldPtr_(NULL)
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
    redirectPatchFieldPtr_(NULL)
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
    redirectPatchFieldPtr_(NULL)
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
        Info<< "constructing patchField from :" << dict << endl;

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
    //dict_.set("value", static_cast<const scalarField&>(*this));
    //os << dict_ << token::END_STATEMENT << nl;
    fixedValueFvPatchField<scalar>::write(os);
    os.writeKeyword("redirectType") << redirectType_ << token::END_STATEMENT
        << nl;
    if (dict_.found("code"))
    {
        os.writeKeyword("code") << string(dict_["code"]) << token::END_STATEMENT
            << nl;
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

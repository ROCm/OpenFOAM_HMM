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


void Foam::codedFixedValueFvPatchScalarField::updateLibrary
(
    bool firstTime
) const
{
    dynamicCode::checkSecurity
    (
        "codedFixedValueFvPatchScalarField::updateLibrary()",
        dict_
    );

    // use codeProperties or in-line
    const bool useCodeProps = !dict_.found("code");

    const dictionary& codeDict =
    (
        useCodeProps
      ? this->dict().subDict(redirectType_)
      : dict_
    );


    autoPtr<dynamicCodeContext> contextPtr;

    // write code into redirectType_ subdir as well
    dynamicCode dynCode(redirectType_);
    const fileName libPath = dynCode.libPath();

    // see if library is loaded
    void* lib = dlLibraryTable::findLibrary(libPath);

    bool reuseLib = false;
    bool waiting = false;

    if (useCodeProps)
    {
        // library may be loaded, but out-of-date
        const codeProperties& codeProps = this->dict();
        if (codeProps.modified())
        {
            codeProps.setUnmodified();

            // Remove instantiation of fvPatchField provided by library
            redirectPatchFieldPtr_.clear();

            contextPtr.reset(new dynamicCodeContext(codeDict));

            // unload code
            if (lib)
            {
                firstTime = false;
                reuseLib = false;
                lib = 0;

                if (!dlLibraryTable::close(libPath, false))
                {
                    FatalIOErrorIn
                    (
                        "codedFixedValueFvPatchScalarField::"
                        "updateLibrary()",
                        contextPtr().dict()
                    )   << "Failed unloading library "
                        << libPath
                        << exit(FatalIOError);
                }
            }
        }
    }


    // library exists (and was not unloaded) - we can leave now
    if (lib)
    {
        return;
    }


    // Remove instantiation of fvPatchField provided by library
    redirectPatchFieldPtr_.clear();

    if (contextPtr.empty())
    {
        contextPtr.reset(new dynamicCodeContext(codeDict));
    }

    // function name serving as version control - based on the SHA1
    const string sentinelName
        = dynCode.codeName() + contextPtr().sha1().str(true);

    // avoid compilation (first time only) by loading an existing library
    if (firstTime && dlLibraryTable::open(libPath, false))
    {
        lib = dlLibraryTable::findLibrary(libPath);

        // verify the loaded version and unload if needed
        if (lib)
        {
            reuseLib = dlSymFound(lib, sentinelName);
            if (!reuseLib)
            {
                lib = 0;
                if (!dlLibraryTable::close(libPath, false))
                {
                    FatalIOErrorIn
                    (
                    "codedFixedValueFvPatchScalarField::updateLibrary()",
                        contextPtr().dict()
                    )   << "Failed unloading library "
                        << libPath
                        << exit(FatalIOError);
                }
            }
        }
    }


    // really do need to create library
    if (!lib)
    {
        if (Pstream::master())
        {
            createLibrary(dynCode, contextPtr());

            if (!dynCode.wmakeLibso())
            {
                FatalIOErrorIn
                (
                    "codedFixedValueFvPatchScalarField::updateLibrary()",
                    contextPtr().dict()
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
                contextPtr().dict()
            )   << "Failed loading library " << libPath
                << exit(FatalIOError);
        }

        lib = dlLibraryTable::findLibrary(libPath);
        if (!lib)
        {
            FatalIOErrorIn
            (
                "codedFixedValueFvPatchScalarField::"
                "updateLibrary()",
                contextPtr().dict()
            )   << "Failed to load library " << libPath
                << exit(FatalIOError);
        }

//#if 0
//        Info<<"check " << libPath << " for " << sentinelName << nl;
//        // paranoid - check that signature function is really there
//        lib = dlLibraryTable::findLibrary(libPath);
//        if (!lib || !dlSymFound(lib, sentinelName))
//        {
//            FatalIOErrorIn
//            (
//                "codedFixedValueFvPatchScalarField::"
//                "updateLibrary()",
//                contextPtr().dict()
//            )   << "Failed to load library with correct signature "
//                << libPath
//                << exit(FatalIOError);
//        }
//#endif
    }
    else if (reuseLib)
    {
        Info<< "Reusing library in " << libPath << nl;
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
    redirectPatchFieldPtr_()
{
    updateLibrary(true);
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
        Info<< "constructing patchField from :" << dict << endl;

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

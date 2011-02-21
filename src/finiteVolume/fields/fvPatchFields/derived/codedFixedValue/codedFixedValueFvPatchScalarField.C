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
#include "codeStreamTools.H"
#include "codeProperties.H"

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


void Foam::codedFixedValueFvPatchScalarField::writeLibrary
(
    const fileName dir,
    const fileName libPath,
    const dictionary& dict
)
{
    Info<< "Creating new library in " << libPath << endl;

    // Write files for new library
    if (Pstream::master())
    {
        fileName templates(Foam::getEnv("OTF_TEMPLATE_DIR"));
        if (!templates.size())
        {
            FatalIOErrorIn
            (
                "codedFixedValueFvPatchScalarField::writeLibrary(..)",
                dict
            )   << "Please set environment variable OTF_TEMPLATE_DIR"
                << " to point to the location of "
                << "fixedValueFvPatchScalarFieldTemplate.C"
                << exit(FatalIOError);
        }


        // Extract sections of code
        string codeInclude = "";
        if (dict.found("codeInclude"))
        {
            codeInclude = codeStreamTools::stripLeading(dict["codeInclude"]);
        }
        string code = codeStreamTools::stripLeading(dict["code"]);

        string codeOptions = "";
        if (dict.found("codeOptions"))
        {
            codeOptions = codeStreamTools::stripLeading(dict["codeOptions"]);
        }


        List<fileAndVars> copyFiles(2);
        copyFiles[0].first() =
            templates/"fixedValueFvPatchScalarFieldTemplate.C";

        copyFiles[0].second().setSize(2);
        copyFiles[0].second()[0] = Pair<string>("OTF_INCLUDES", codeInclude);
        copyFiles[0].second()[1] = Pair<string>("OTF_UPDATECOEFFS", code);

        copyFiles[1].first() =
            templates/"fixedValueFvPatchScalarFieldTemplate.H";



        List<fileAndContent> filesContents(2);
        // Write Make/files
        filesContents[0].first() = "Make/files";
        filesContents[0].second() =
            "fixedValueFvPatchScalarFieldTemplate.C \n\n"
            "LIB = $(FOAM_USER_LIBBIN)/lib" + redirectType_;
        // Write Make/options
        filesContents[1].first() = "Make/options";
        filesContents[1].second() =
            "EXE_INC = -g\\\n    -I$(LIB_SRC)/finiteVolume/lnInclude\\\n"
          + codeOptions
          + "\n\nLIB_LIBS = ";

        codeStreamTools writer(redirectType_, copyFiles, filesContents);
        if (!writer.copyFilesContents(dir))
        {
            FatalIOErrorIn
            (
                "codedFixedValueFvPatchScalarField::writeLibrary(..)",
                dict
            )   << "Failed writing " << endl
                << copyFiles << endl
                << filesContents
                << exit(FatalIOError);
        }
    }
}


void Foam::codedFixedValueFvPatchScalarField::updateLibrary()
{
    if (isAdministrator())
    {
        FatalIOErrorIn
        (
            "codedFixedValueFvPatchScalarField::updateLibrary()",
            dict_
        )   << "This code should not be executed by someone with administrator"
            << " rights due to security reasons." << endl
            << "(it writes a shared library which then gets loaded "
            << "using dlopen)"
            << exit(FatalIOError);
    }

    const fileName dir =
        db().time().constantPath()/"codedFixedValue"/redirectType_;
    //Info<< "dir:" << dir << endl;

    const fileName libPath
    (
        Foam::getEnv("FOAM_USER_LIBBIN")
      / "lib"
      + redirectType_
      + ".so"
    );
    //Info<< "libPath:" << libPath << endl;

    void* lib = dlLibraryTable::findLibrary(libPath);

    if (dict_.found("code"))
    {
        if (!lib)
        {
            writeLibrary(dir, libPath, dict_);
        }
    }
    else
    {
        const codeProperties& onTheFlyDict = dict();

        if (onTheFlyDict.modified())
        {
            onTheFlyDict.setUnmodified();

            // Remove instantiation of fvPatchField provided by library
            redirectPatchFieldPtr_.clear();
            // Unload library
            if (lib)
            {
                if (!dlLibraryTable::close(libPath))
                {
                    FatalIOErrorIn
                    (
                        "codedFixedValueFvPatchScalarField::updateLibrary(..)",
                        onTheFlyDict
                    )   << "Failed unloading library " << libPath
                        << exit(FatalIOError);
                }
                lib = NULL;
            }

            const dictionary& codeDict = onTheFlyDict.subDict(redirectType_);
            writeLibrary(dir, libPath, codeDict);
        }
    }

    if (!lib)
    {
        if (Pstream::master())
        {
            Foam::string wmakeCmd("wmake libso " + dir);
            Info<< "Invoking " << wmakeCmd << endl;
            if (Foam::system(wmakeCmd))
            {
                FatalIOErrorIn
                (
                    "codedFixedValueFvPatchScalarField::updateLibrary()",
                    dict_
                )   << "Failed " << wmakeCmd << exit(FatalIOError);
            }
        }

        bool dummy = true;
        reduce(dummy, orOp<bool>());

        if (!dlLibraryTable::open(libPath))
        {
            FatalIOErrorIn
            (
                "codedFixedValueFvPatchScalarField::updateLibrary()",
                dict_
            )   << "Failed loading library " << libPath << exit(FatalIOError);
        }
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

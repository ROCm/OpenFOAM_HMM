/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "yPlusRAS.H"
#include "volFields.H"

#include "incompressible/RAS/RASModel/RASModel.H"
#include "nutkWallFunction/nutkWallFunctionFvPatchScalarField.H"

#include "compressible/RAS/RASModel/RASModel.H"
#include "mutkWallFunction/mutkWallFunctionFvPatchScalarField.H"

#include "wallDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::yPlusRAS, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::yPlusRAS::makeFile()
{
    // Create the output file if not already created
    if (outputFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating output file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName outputDir;
            word startTimeName =
                obr_.time().timeName(obr_.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                outputDir =
                    obr_.time().path()/".."/name_/startTimeName;
            }
            else
            {
                outputDir = obr_.time().path()/name_/startTimeName;
            }

            // Create directory if does not exist
            mkDir(outputDir);

            // Open new file at start up
            outputFilePtr_.reset(new OFstream(outputDir/(type() + ".dat")));

            // Add headers to output data
            outputFilePtr_() << "# y+ (RAS)" << nl
                << "# time " << token::TAB << "patch" << token::TAB
                << "min" << token::TAB << "max" << token::TAB << "average"
                << endl;
        }
    }
}


void Foam::yPlusRAS::calcIncompressibleYPlus
(
    const fvMesh& mesh,
    volScalarField& yPlus
)
{
    typedef incompressible::RASModels::nutkWallFunctionFvPatchScalarField
        wallFunctionPatchField;

    const incompressible::RASModel& model =
        mesh.lookupObject<incompressible::RASModel>("RASProperties");

    const volScalarField nut(model.nut());
    const volScalarField::GeometricBoundaryField& nutPatches =
        nut.boundaryField();

    bool foundPatch = false;
    forAll(nutPatches, patchI)
    {
        if (isA<wallFunctionPatchField>(nutPatches[patchI]))
        {
            foundPatch = true;

            const wallFunctionPatchField& nutPw =
                dynamic_cast<const wallFunctionPatchField&>(nutPatches[patchI]);

            yPlus.boundaryField()[patchI] = nutPw.yPlus();
            const scalarField& Yp = yPlus.boundaryField()[patchI];

            scalar minYp = min(Yp);
            scalar maxYp = max(Yp);
            scalar avgYp = average(Yp);

            if (log_)
            {
                Info<< "    patch " << nutPw.patch().name()
                    << " y+ : min = " << minYp << ", max = " << maxYp
                    << ", average = " << avgYp << nl;
            }

            if (Pstream::master())
            {
                outputFilePtr_() << obr_.time().value() << token::TAB
                    << nutPw.patch().name() << token::TAB
                    << minYp << token::TAB << maxYp << token::TAB
                    << avgYp << endl;
            }
        }
    }

    if (log_ && !foundPatch)
    {
        Info<< "    no " << wallFunctionPatchField::typeName << " patches"
            << endl;
    }
}


void Foam::yPlusRAS::calcCompressibleYPlus
(
    const fvMesh& mesh,
    volScalarField& yPlus
)
{
    typedef compressible::RASModels::mutkWallFunctionFvPatchScalarField
        wallFunctionPatchField;

    const compressible::RASModel& model =
        mesh.lookupObject<compressible::RASModel>("RASProperties");

    const volScalarField mut(model.mut());
    const volScalarField::GeometricBoundaryField& mutPatches =
        mut.boundaryField();

    bool foundPatch = false;
    forAll(mutPatches, patchI)
    {
        if (isA<wallFunctionPatchField>(mutPatches[patchI]))
        {
            foundPatch = true;

            const wallFunctionPatchField& mutPw =
                dynamic_cast<const wallFunctionPatchField&>(mutPatches[patchI]);

            yPlus.boundaryField()[patchI] = mutPw.yPlus();
            const scalarField& Yp = yPlus.boundaryField()[patchI];

            scalar minYp = min(Yp);
            scalar maxYp = max(Yp);
            scalar avgYp = average(Yp);

            if (log_)
            {
                Info<< "    patch " << mutPw.patch().name()
                    << " y+ : min = " << minYp << ", max = " << maxYp
                    << ", average = " << avgYp << nl;
            }

            if (Pstream::master())
            {
                outputFilePtr_() << obr_.time().value() << token::TAB
                    << mutPw.patch().name() << token::TAB
                    << minYp << token::TAB << maxYp << token::TAB
                    << avgYp << endl;
            }
        }
    }

    if (log_ && !foundPatch)
    {
        Info<< "    no " << wallFunctionPatchField::typeName << " patches"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::yPlusRAS::yPlusRAS
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    log_(false),
    phiName_("phi"),
    outputFilePtr_(NULL)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "yPlusRAS::yPlusRAS"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }

    makeFile();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::yPlusRAS::~yPlusRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::yPlusRAS::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);
        phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
    }
}


void Foam::yPlusRAS::execute()
{
    // Do nothing - only valid on write
}


void Foam::yPlusRAS::end()
{
    // Do nothing - only valid on write
}


void Foam::yPlusRAS::write()
{
    if (active_)
    {
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(phiName_);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField yPlusRAS
        (
            IOobject
            (
                "yPlusRAS",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar("0", dimless, 0.0)
        );

        if (log_)
        {
            Info<< type() << " output:" << nl;
        }

        if (phi.dimensions() == dimMass/dimTime)
        {
            calcCompressibleYPlus(mesh, yPlusRAS);
        }
        else
        {
            calcIncompressibleYPlus(mesh, yPlusRAS);
        }

        if (log_)
        {
            Info<< endl;
        }

        yPlusRAS.write();
    }
}


// ************************************************************************* //

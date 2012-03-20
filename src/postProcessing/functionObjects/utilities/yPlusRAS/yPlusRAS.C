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

void Foam::yPlusRAS::calcIncompressibleYPlus
(
    const fvMesh& mesh,
    volScalarField& yPlus
) const
{
    typedef incompressible::RASModels::nutkWallFunctionFvPatchScalarField
        wallFunctionPatchField;

    const incompressible::RASModel& model =
        mesh.lookupObject<incompressible::RASModel>("RASProperties");

    const volScalarField nut(model.nut());
    const volScalarField::GeometricBoundaryField& nutPatches =
        nut.boundaryField();

    Info<< type() << " output:" << nl;

    bool foundNutPatch = false;
    forAll(nutPatches, patchi)
    {
        if (isA<wallFunctionPatchField>(nutPatches[patchi]))
        {
            foundNutPatch = true;

            const wallFunctionPatchField& nutPw =
                dynamic_cast<const wallFunctionPatchField&>(nutPatches[patchi]);

            yPlus.boundaryField()[patchi] = nutPw.yPlus();
            const scalarField& Yp = yPlus.boundaryField()[patchi];

            Info<< "    patch " << nutPw.patch().name()
                << " y+ : min = " << min(Yp) << ", max = " << max(Yp)
                << ", average = " << average(Yp) << nl;
        }
    }

    if (!foundNutPatch)
    {
        Info<< "    no " << wallFunctionPatchField::typeName << " patches"
            << endl;
    }
}


void Foam::yPlusRAS::calcCompressibleYPlus
(
    const fvMesh& mesh,
    volScalarField& yPlus
) const
{
    typedef compressible::RASModels::mutkWallFunctionFvPatchScalarField
        wallFunctionPatchField;

    const compressible::RASModel& model =
        mesh.lookupObject<compressible::RASModel>("RASProperties");

    const volScalarField mut(model.mut());
    const volScalarField::GeometricBoundaryField& mutPatches =
        mut.boundaryField();

    Info<< type() << " output:" << nl;

    bool foundMutPatch = false;
    forAll(mutPatches, patchi)
    {
        if (isA<wallFunctionPatchField>(mutPatches[patchi]))
        {
            foundMutPatch = true;

            const wallFunctionPatchField& mutPw =
                dynamic_cast<const wallFunctionPatchField&>(mutPatches[patchi]);

            yPlus.boundaryField()[patchi] = mutPw.yPlus();
            const scalarField& Yp = yPlus.boundaryField()[patchi];

            Info<< "    patch " << mutPw.patch().name()
                << " y+ : min = " << min(Yp) << ", max = " << max(Yp)
                << ", average = " << average(Yp) << nl;
        }
    }

    if (!foundMutPatch)
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
    phiName_("phi")
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::yPlusRAS::~yPlusRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::yPlusRAS::read(const dictionary& dict)
{
    if (active_)
    {
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

        if (phi.dimensions() == dimMass/dimTime)
        {
            calcCompressibleYPlus(mesh, yPlusRAS);
        }
        else
        {
            calcIncompressibleYPlus(mesh, yPlusRAS);
        }

        Info<< endl;

        yPlusRAS.write();
    }
}


// ************************************************************************* //

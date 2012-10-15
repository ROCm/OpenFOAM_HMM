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

#include "wallShearStress.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::wallShearStress, 0);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::wallShearStress::writeFileHeader(const label i)
{
    // Add headers to output data
    file() << "# Wall shear stress" << nl
        << "# time " << token::TAB << "patch" << token::TAB
        << "min" << token::TAB << "max" << endl;
}


void Foam::wallShearStress::calcShearStress
(
    const fvMesh& mesh,
    const volSymmTensorField& Reff,
    volVectorField& shearStress
)
{
    forAll(shearStress.boundaryField(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        if (isA<wallPolyPatch>(pp))
        {
            vectorField& ssp = shearStress.boundaryField()[patchI];
            const vectorField& Sfp = mesh.Sf().boundaryField()[patchI];
            const scalarField& magSfp = mesh.magSf().boundaryField()[patchI];
            const symmTensorField& Reffp = Reff.boundaryField()[patchI];

            ssp = (-Sfp/magSfp) & Reffp;

            vector minSsp = min(ssp);
            vector maxSsp = max(ssp);

            file() << mesh.time().timeName() << token::TAB
                << pp.name() << token::TAB << minSsp
                << token::TAB << maxSsp << endl;

            if (log_)
            {
                Info<< "    min/max(" << pp.name() << ") = "
                    << minSsp << ", " << maxSsp << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallShearStress::wallShearStress
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    obr_(obr),
    active_(true),
    log_(false),
    phiName_("phi")
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "wallShearStress::wallShearStress"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallShearStress::~wallShearStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallShearStress::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);
        phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
    }
}


void Foam::wallShearStress::execute()
{
    // Do nothing - only valid on write
}


void Foam::wallShearStress::end()
{
    // Do nothing - only valid on write
}


void Foam::wallShearStress::write()
{
    typedef compressible::turbulenceModel cmpModel;
    typedef incompressible::turbulenceModel icoModel;

    if (active_)
    {
        functionObjectFile::write();

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volVectorField wallShearStress
        (
            IOobject
            (
                "wallShearStress",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedVector("0", sqr(dimLength)/sqr(dimTime), vector::zero)
        );

        if (log_)
        {
            Info<< type() << " output:" << nl;
        }


        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(phiName_);

        tmp<volSymmTensorField> Reff;
        if (phi.dimensions() == dimMass/dimTime)
        {
            if (!mesh.foundObject<cmpModel>("turbulenceModel"))
            {
                FatalErrorIn("void Foam::wallShearStress::write()")
                    << "Unable to find compressible turbulence model in the "
                    << "database" << exit(FatalError);
            }

            const cmpModel& model =
                mesh.lookupObject<cmpModel>("turbulenceModel");

            Reff = model.devRhoReff();
        }
        else
        {
            if (!mesh.foundObject<icoModel>("turbulenceModel"))
            {
                FatalErrorIn("void Foam::wallShearStress::write()")
                    << "Unable to find incompressible turbulence model in the "
                    << "database" << exit(FatalError);
            }

            const icoModel& model =
                mesh.lookupObject<icoModel>("turbulenceModel");

            Reff = model.devReff();
        }

    
        calcShearStress(mesh, Reff(), wallShearStress);

        if (log_)
        {
            Info<< endl;
        }

        wallShearStress.write();
    }
}


// ************************************************************************* //

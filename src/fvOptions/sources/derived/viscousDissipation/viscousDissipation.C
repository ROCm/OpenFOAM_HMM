/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "viscousDissipation.H"
#include "fvMatrices.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(viscousDissipation, 0);

    addToRunTimeSelectionTable
    (
        option,
        viscousDissipation,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fv::viscousDissipation::rho() const
{
    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                "trho",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            rho_
        )
    );

    if (rho_.value() > 0)
    {
        return trho;
    }
    else if (rhoName_ != "none")
    {
        trho.ref() = mesh_.lookupObject<volScalarField>(rhoName_);
        return trho;
    }

    FatalErrorInFunction
        << "Neither rhoName nor rho are specified."
        << exit(FatalError);

    return tmp<volScalarField>();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::viscousDissipation::viscousDissipation
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(sourceName, modelType, dict, mesh),
    UName_(coeffs_.lookupOrDefault<word>("U", "U")),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "none")),
    rho_
    (
        coeffs_.lookupOrDefault
        (
            "rhoInf",
            dimensionedScalar("rho", dimDensity, 0)
        )
    )
{
    const basicThermo* thermoPtr =
        mesh_.lookupObjectPtr<basicThermo>(basicThermo::dictName);

    if (thermoPtr)
    {
        fieldNames_.setSize(1, thermoPtr->he().name());
    }

    if (fieldNames_.empty())
    {
        coeffs_.lookup("fields") >> fieldNames_;
    }

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    applied_.setSize(fieldNames_.size(), false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField> Foam::fv::viscousDissipation::
devRhoReff() const
{
    // Incompressible
    {
        typedef incompressible::turbulenceModel turbType;

        const turbType* turbPtr =
            mesh_.lookupObjectPtr<turbType>(turbulenceModel::propertiesName);

        if (turbPtr)
        {
            return tmp<volSymmTensorField>(rho()*turbPtr->devRhoReff());
        }
    }

    // Compressible
    {
        typedef compressible::turbulenceModel turbType;

        const turbType* turbPtr =
            mesh_.lookupObjectPtr<turbType>(turbulenceModel::propertiesName);

        if (turbPtr)
        {
            return tmp<volSymmTensorField>(turbPtr->devRhoReff());
        }
    }

    FatalErrorInFunction
        << " The turbulence model is not found in the database."
        << exit(FatalError);

    return tmp<volSymmTensorField>();
}


void Foam::fv::viscousDissipation::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    typedef typename outerProduct<vector, vector>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    word gradUName("grad(" + UName_ + ')');

    tmp<GradFieldType> tgradU
    (
        new GradFieldType
        (
            IOobject
            (
                "gradU",
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor("zero", inv(dimTime) , tensor::zero)
        )
    );

    // Cached?
    const GradFieldType* gradUPtr =
        mesh_.lookupObjectPtr<GradFieldType>(gradUName);

    if (gradUPtr)
    {
        tgradU.ref() = *gradUPtr;
    }
    else
    {
        const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
        tgradU.ref() = fvc::grad(U);
    }

    const volScalarField D("D", devRhoReff() && tgradU.ref());

    eqn -= D;
}


// ************************************************************************* //

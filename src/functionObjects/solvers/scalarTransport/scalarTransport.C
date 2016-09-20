/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "scalarTransport.H"
#include "surfaceFields.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(scalarTransport, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        scalarTransport,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::volScalarField& Foam::scalarTransport::transportedField()
{
    if (!mesh_.foundObject<volScalarField>(fieldName_))
    {
        volScalarField* fldPtr = new volScalarField
        (
            IOobject
            (
                fieldName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0),
            boundaryTypes()
        );
        fldPtr->store();
    }

    return const_cast<volScalarField&>
    (
        mesh_.lookupObject<volScalarField>(fieldName_)
    );
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::scalarTransport::D
(
    const surfaceScalarField& phi
) const
{
    typedef incompressible::turbulenceModel icoModel;
    typedef compressible::turbulenceModel cmpModel;

    word Dname("D" + s_.name());

    if (constantD_)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    Dname,
                    mesh_.time().timeName(),
                    mesh_.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(Dname, phi.dimensions()/dimLength, D_)
            )
        );
    }
    else if (mesh_.foundObject<icoModel>(turbulenceModel::propertiesName))
    {
        const icoModel& model = mesh_.lookupObject<icoModel>
        (
            turbulenceModel::propertiesName
        );

        return model.nuEff();
    }
    else if (mesh_.foundObject<cmpModel>(turbulenceModel::propertiesName))
    {
        const cmpModel& model = mesh_.lookupObject<cmpModel>
        (
            turbulenceModel::propertiesName
        );

        return model.muEff();
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    Dname,
                    mesh_.time().timeName(),
                    mesh_.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(Dname, phi.dimensions()/dimLength, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::scalarTransport::scalarTransport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_(dict.lookupOrDefault<word>("field", "s")),
    phiName_("phi"),
    rhoName_("rho"),
    D_(0),
    constantD_(false),
    nCorr_(0),
    resetOnStartUp_(false),
    schemesField_("unknown-schemesField"),
    fvOptions_(mesh_)
{
    read(dict);

    // Force creation of transported field so any BCs using it can
    // look it up
    volScalarField& s = transportedField();

    if (resetOnStartUp_)
    {
        s == dimensionedScalar("zero", dimless, 0.0);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::scalarTransport::~scalarTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::scalarTransport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("phi", phiName_);
    dict.readIfPresent("rho", rhoName_);
    schemesField_ = dict.lookupOrDefault("schemesField", fieldName_);

    constantD_ = false;
    if (dict.readIfPresent("D", D_))
    {
        constantD_ = true;
    }

    dict.readIfPresent("nCorr", nCorr_);
    dict.readIfPresent("resetOnStartUp", resetOnStartUp_);

    if (dict.found("fvOptions"))
    {
        fvOptions_.reset(dict.subDict("fvOptions"));
    }

    return true;
}


bool Foam::functionObjects::scalarTransport::execute()
{
    Log << type() << " write:" << endl;

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    volScalarField& s = transportedField();

    // Calculate the diffusivity
    volScalarField D(this->D(phi));

    word divScheme("div(phi," + schemesField_ + ")");
    word laplacianScheme("laplacian(" + D.name() + "," + schemesField_ + ")");

    // Set under-relaxation coeff
    scalar relaxCoeff = 0.0;
    if (mesh_.relaxEquation(schemesField_))
    {
        relaxCoeff = mesh_.equationRelaxationFactor(schemesField_);
    }

    if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        for (label i = 0; i <= nCorr_; i++)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(rho, s)
              + fvm::div(phi, s, divScheme)
              - fvm::laplacian(D, s, laplacianScheme)
             ==
                fvOptions_(rho, s)
            );

            sEqn.relax(relaxCoeff);

            fvOptions_.constrain(sEqn);

            sEqn.solve(mesh_.solverDict(schemesField_));
        }
    }
    else if (phi.dimensions() == dimVolume/dimTime)
    {
        for (label i = 0; i <= nCorr_; i++)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(s)
              + fvm::div(phi, s, divScheme)
              - fvm::laplacian(D, s, laplacianScheme)
             ==
                fvOptions_(s)
            );

            sEqn.relax(relaxCoeff);

            fvOptions_.constrain(sEqn);

            sEqn.solve(mesh_.solverDict(schemesField_));
        }
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << endl;
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::scalarTransport::write()
{
    return true;
}


// ************************************************************************* //

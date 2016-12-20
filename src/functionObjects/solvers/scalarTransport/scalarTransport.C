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
#include "CMULES.H"
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

Foam::volScalarField& Foam::functionObjects::scalarTransport::transportedField()
{
    if (!foundObject<volScalarField>(fieldName_))
    {
        tmp<volScalarField> tfldPtr
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
        store(fieldName_, tfldPtr);

        if (phaseName_ != "none")
        {
            mesh_.setFluxRequired(fieldName_);
        }
    }

    return const_cast<volScalarField&>
    (
        lookupObject<volScalarField>(fieldName_)
    );
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::scalarTransport::D
(
    const volScalarField& s,
    const surfaceScalarField& phi
) const
{
    typedef incompressible::turbulenceModel icoModel;
    typedef compressible::turbulenceModel cmpModel;

    word Dname("D" + s.name());

    if (constantD_)
    {
        tmp<volScalarField> tD
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

        return tD;
    }
    else if (nutName_ != "none")
    {
        const volScalarField& nutMean =
            mesh_.lookupObject<volScalarField>(nutName_);

        return tmp<volScalarField>
        (
            new volScalarField(Dname, nutMean)
        );
    }
    else if (foundObject<icoModel>(turbulenceModel::propertiesName))
    {
        const icoModel& model = lookupObject<icoModel>
        (
            turbulenceModel::propertiesName
        );

        return tmp<volScalarField>
        (
             new volScalarField(Dname, model.nuEff())
        );
    }
    else if (foundObject<cmpModel>(turbulenceModel::propertiesName))
    {
        const cmpModel& model = lookupObject<cmpModel>
        (
            turbulenceModel::propertiesName
        );

        return tmp<volScalarField>
        (
             new volScalarField(Dname, model.muEff())
        );
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
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    UPhiName_(dict.lookupOrDefault<word>("UPhi", "none")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    nutName_(dict.lookupOrDefault<word>("nut", "none")),
    phaseName_(dict.lookupOrDefault<word>("phase", "none")),
    phasePhiCompressedName_
    (
        dict.lookupOrDefault<word>("phasePhiCompressed", "alphaPhiUn")
    ),
    D_(0),
    constantD_(false),
    nCorr_(0),
    resetOnStartUp_(false),
    schemesField_("unknown-schemesField"),
    fvOptions_(mesh_),
    bounded01_(dict.lookupOrDefault<bool>("bounded01", true))
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
    dict.readIfPresent("UPhi", UPhiName_);
    dict.readIfPresent("nut", nutName_);
    dict.readIfPresent("phase", phaseName_);
    dict.readIfPresent("bounded01", bounded01_);

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

    tmp<surfaceScalarField> tPhi
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("tPhi", dimMass/dimTime, 0.0)
        )
    );
    surfaceScalarField& phi = tPhi.ref();

    const dimensionSet dim
    (
        mesh_.lookupObject<surfaceScalarField>(phiName_).dimensions()
    );

    bool compressible = true;
    if (dim == dimVolume/dimTime)
    {
        compressible = false;
        phi.dimensions().reset(dimVolume/dimTime);
    }

    // Obtain phi from phiName or constructed from UPhiName
    if (phiName_ != "none")
    {
        phi = const_cast<surfaceScalarField&>
        (
            mesh_.lookupObject<surfaceScalarField>(phiName_)
        );
    }
    else if (UPhiName_ != "none")
    {
        const volVectorField& Uphi =
            mesh_.lookupObject<volVectorField>(UPhiName_);

        if (!compressible)
        {
            phi = fvc::interpolate(Uphi) & mesh_.Sf();
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName_);

            phi = fvc::interpolate(rho*Uphi) & mesh_.Sf();
        }
    }

    volScalarField& s = transportedField();

    // Calculate the diffusivity
    volScalarField D(this->D(s, phi));

    word divScheme("div(phi," + schemesField_ + ")");
    word laplacianScheme("laplacian(" + D.name() + "," + schemesField_ + ")");

    // Set under-relaxation coeff
    scalar relaxCoeff = 0.0;
    if (mesh_.relaxEquation(schemesField_))
    {
        relaxCoeff = mesh_.equationRelaxationFactor(schemesField_);
    }

    // Two phase scalar transport
    if (phaseName_ != "none")
    {
        const volScalarField& alpha =
            mesh_.lookupObject<volScalarField>(phaseName_);

        const surfaceScalarField& limitedPhiAlpa =
            mesh_.lookupObject<surfaceScalarField>(phasePhiCompressedName_);

        D *= pos(alpha - 0.99);
/*
        surfaceScalarField phic(2.0*mag(phi/mesh_.magSf()));

        const volVectorField gradAlpha(fvc::grad(alpha, "nHat"));

        surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

        dimensionedScalar deltaN
        (
            "deltaN", 1e-8/pow(average(mesh_.V()), 1.0/3.0)
        );

        surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN));

        surfaceScalarField nHat(nHatfv & mesh_.Sf());

        surfaceScalarField phir(phic*nHat);

        surfaceScalarField limitedPhiAlpa
        (
            fvc::flux
            (
                phi,
                alpha,
                "div(phi,s)"
            )
          + fvc::flux
            (
               -fvc::flux(-phir, (1-alpha), "div(phirb,s)"),
                alpha,
                "div(phirb,s)"
            )
        );
*/
        // Reset D dimensions consistent with limitedPhiAlpa
        D.dimensions().reset(limitedPhiAlpa.dimensions()/dimLength);

        // Solve
        tmp<surfaceScalarField> tTPhiUD;
        for (label i = 0; i <= nCorr_; i++)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(s)
              + fvm::div(limitedPhiAlpa, s, divScheme)
              - fvm::laplacian(D, s, laplacianScheme)
              ==
                alpha*fvOptions_(s)
            );

            sEqn.relax(relaxCoeff);
            fvOptions_.constrain(sEqn);
            sEqn.solve(mesh_.solverDict(schemesField_));

            tTPhiUD = sEqn.flux();
        }

        if (bounded01_)
        {
            MULES::explicitSolve(s, phi, tTPhiUD.ref(), 1, 0);
        }
    }
    else if (compressible)
    {
        const volScalarField& rho = lookupObject<volScalarField>(rhoName_);

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
    else if (!compressible)
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

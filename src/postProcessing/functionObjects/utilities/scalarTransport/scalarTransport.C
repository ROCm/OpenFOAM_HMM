/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
#include "dictionary.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fvScalarMatrix.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(scalarTransport, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::scalarTransport::boundaryTypes() const
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    wordList bTypes(U.boundaryField().size());

    forAll(bTypes, patchI)
    {
        const fvPatchField<vector>& pf = U.boundaryField()[patchI];
        if (isA<fixedValueFvPatchVectorField>(pf))
        {
            bTypes[patchI] = fixedValueFvPatchScalarField::typeName;
        }
        else
        {
            bTypes[patchI] = zeroGradientFvPatchScalarField::typeName;
        }
    }

    return bTypes;
}


Foam::volScalarField& Foam::scalarTransport::transportedField()
{
    if (!mesh_.foundObject<volScalarField>(name()))
    {
        volScalarField* fldPtr = new volScalarField
        (
            IOobject
            (
                name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
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
        mesh_.lookupObject<volScalarField>(name())
    );
}


Foam::tmp<Foam::volScalarField> Foam::scalarTransport::DT
(
    const surfaceScalarField& phi
) const
{
    typedef incompressible::turbulenceModel icoModel;
    typedef compressible::turbulenceModel cmpModel;

    if (userDT_)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "DT",
                    mesh_.time().timeName(),
                    mesh_.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("DT", phi.dimensions()/dimLength, DT_)
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
                    "DT",
                    mesh_.time().timeName(),
                    mesh_.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("DT", phi.dimensions()/dimLength, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarTransport::scalarTransport
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    mesh_(refCast<const fvMesh>(obr)),
    active_(true),
    phiName_("phi"),
    UName_("U"),
    rhoName_("rho"),
    DT_(0.0),
    userDT_(false),
    resetOnStartUp_(false),
    nCorr_(0),
    autoSchemes_(false),
    fvOptions_(mesh_),
    log_(true)
{
    read(dict);

    // Force creation of transported field so any bcs using it can look it
    // up
    volScalarField& T = transportedField();

    if (resetOnStartUp_)
    {
        T == dimensionedScalar("zero", dimless, 0.0);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::scalarTransport::~scalarTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::scalarTransport::read(const dictionary& dict)
{
    if (active_)
    {
        log_.readIfPresent("log", dict);

        if (log_) Info<< type() << " " << name_ << " output:" << nl;

        dict.readIfPresent("phiName", phiName_);
        dict.readIfPresent("UName", UName_);
        dict.readIfPresent("rhoName", rhoName_);

        userDT_ = false;
        if (dict.readIfPresent("DT", DT_))
        {
            userDT_ = true;
        }

        dict.readIfPresent("nCorr", nCorr_);
        dict.lookup("resetOnStartUp") >> resetOnStartUp_;
        dict.lookup("autoSchemes") >> autoSchemes_;

        fvOptions_.reset(dict.subDict("fvOptions"));
    }
}


void Foam::scalarTransport::execute()
{
    if (active_)
    {
        if (log_) Info<< type() << " " << name_ << " output:" << nl;

        const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>(phiName_);

        volScalarField& T = transportedField();

        // calculate the diffusivity
        volScalarField DT(this->DT(phi));

        // set schemes
        word schemeVar = T.name();
        if (autoSchemes_)
        {
            schemeVar = UName_;
        }

        word divScheme("div(phi," + schemeVar + ")");
        word laplacianScheme("laplacian(" + DT.name() + "," + schemeVar + ")");

        // set under-relaxation coeff
        scalar relaxCoeff = 0.0;
        if (mesh_.relaxEquation(schemeVar))
        {
            relaxCoeff = mesh_.equationRelaxationFactor(schemeVar);
        }

        if (phi.dimensions() == dimMass/dimTime)
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName_);

            // solve
            for (label i = 0; i <= nCorr_; i++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(rho, T)
                  + fvm::div(phi, T, divScheme)
                  - fvm::laplacian(DT, T, laplacianScheme)
                 ==
                    fvOptions_(rho, T)
                );

                TEqn.relax(relaxCoeff);

                fvOptions_.constrain(TEqn);

                TEqn.solve(mesh_.solverDict(schemeVar));
            }
        }
        else if (phi.dimensions() == dimVolume/dimTime)
        {
            // solve
            for (label i = 0; i <= nCorr_; i++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  + fvm::div(phi, T, divScheme)
                  - fvm::laplacian(DT, T, laplacianScheme)
                 ==
                    fvOptions_(T)
                );

                TEqn.relax(relaxCoeff);

                fvOptions_.constrain(TEqn);

                TEqn.solve(mesh_.solverDict(schemeVar));
            }
        }
        else
        {
            FatalErrorInFunction
                << "Incompatible dimensions for phi: " << phi.dimensions() << nl
                << "Dimensions should be " << dimMass/dimTime << " or "
                << dimVolume/dimTime << endl;
        }

        if (log_) Info<< endl;
    }
}


void Foam::scalarTransport::end()
{
    // Do nothing
}


void Foam::scalarTransport::timeSet()
{
    // Do nothing
}


void Foam::scalarTransport::write()
{
    // Do nothing
}


// ************************************************************************* //

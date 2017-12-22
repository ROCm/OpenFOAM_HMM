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

#include "ReynoldsAnalogy.H"
#include "fluidThermo.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{
    defineTypeNameAndDebug(ReynoldsAnalogy, 0);
    addToRunTimeSelectionTable
    (
        heatTransferCoeffModel,
        ReynoldsAnalogy,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::heatTransferCoeffModels::ReynoldsAnalogy::rho(const label patchi) const
{
    if (rhoName_ == "rhoInf")
    {
        const label n = mesh_.boundary()[patchi].size();
        return tmp<Field<scalar>>(new Field<scalar>(n, rhoRef_));
    }
    else if (mesh_.foundObject<volScalarField>(rhoName_, false))
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);
        return rho.boundaryField()[patchi];
    }
    else
    {
        FatalErrorInFunction
            << "Unable to set rho for patch " << patchi
            << exit(FatalError);
    }

    return tmp<Field<scalar>>(nullptr);
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::heatTransferCoeffModels::ReynoldsAnalogy::Cp(const label patchi) const
{
    if (CpName_ == "CpInf")
    {
        const label n = mesh_.boundary()[patchi].size();
        return tmp<Field<scalar>>(new Field<scalar>(n, CpRef_));
    }
    else if (mesh_.foundObject<fluidThermo>(fluidThermo::typeName))
    {
        const fluidThermo& thermo =
            mesh_.lookupObject<fluidThermo>(fluidThermo::typeName);

        const scalarField& pp = thermo.p().boundaryField()[patchi];
        const scalarField& Tp = thermo.T().boundaryField()[patchi];

        return thermo.Cp(pp, Tp, patchi);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to set Cp for patch " << patchi
            << exit(FatalError);
    }

    return tmp<Field<scalar>>(nullptr);
}


Foam::tmp<Foam::volSymmTensorField>
Foam::heatTransferCoeffModels::ReynoldsAnalogy::devReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            mesh_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff()/turb.rho();
    }
    else if (mesh_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const incompressible::turbulenceModel& turb =
            mesh_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return turb.devReff();
    }
    else if (mesh_.foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo& thermo =
            mesh_.lookupObject<fluidThermo>(fluidThermo::dictName);

        const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

        return -thermo.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (mesh_.foundObject<transportModel>("transportProperties"))
    {
        const transportModel& laminarT =
            mesh_.lookupObject<transportModel>("transportProperties");

        const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

        return -laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (mesh_.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
            mesh_.lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu(transportProperties.lookup("nu"));

        const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

        return -nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


Foam::tmp<Foam::FieldField<Foam::Field, Foam::scalar>>
Foam::heatTransferCoeffModels::ReynoldsAnalogy::Cf() const
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    const volVectorField::Boundary& Ubf = U.boundaryField();

    tmp<FieldField<Field, scalar>> tCf
    (
        new FieldField<Field, scalar>(Ubf.size())
    );

    FieldField<Field, scalar>& Cf = tCf.ref();

    forAll(Cf, patchi)
    {
        Cf.set(patchi, new Field<scalar>(Ubf[patchi].size(), 0));
    }

    const volSymmTensorField R(devReff());
    const volSymmTensorField::Boundary& Rbf = R.boundaryField();

    for (label patchi : patchSet_)
    {
        const fvPatchVectorField& Up = Ubf[patchi];

        const symmTensorField& Rp = Rbf[patchi];

        const vectorField nHat(Up.patch().nf());

        const scalarField tauByRhop(mag(nHat & Rp));

        Cf[patchi] = 2*tauByRhop/magSqr(URef_);
    }

    return tCf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferCoeffModels::ReynoldsAnalogy::ReynoldsAnalogy
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& TName
)
:
    heatTransferCoeffModel(dict, mesh, TName),
    UName_("U"),
    URef_(vector::zero),
    rhoName_("rho"),
    rhoRef_(0.0),
    CpName_("Cp"),
    CpRef_(0.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::heatTransferCoeffModels::ReynoldsAnalogy::read
(
    const dictionary& dict
)
{
    if (heatTransferCoeffModel::read(dict))
    {
        dict.lookup("UInf") >> URef_;

        dict.readIfPresent("Cp", CpName_);
        if (CpName_ == "CpInf")
        {
            dict.lookup("CpInf") >> CpRef_;
        }

        dict.readIfPresent("rho", rhoName_);
        if (rhoName_ == "rhoInf")
        {
            dict.lookup("rhoInf") >> rhoRef_;
        }

        return true;
    }

    return false;
}


void Foam::heatTransferCoeffModels::ReynoldsAnalogy::htc(volScalarField& htc)
{
    const FieldField<Field, scalar> CfBf(Cf());
    const scalar magU = mag(URef_);

    volScalarField::Boundary& htcBf = htc.boundaryFieldRef();
    forAllConstIters(patchSet_, iter)
    {
        label patchi = iter.key();
        const scalarField rhop(rho(patchi));
        const scalarField Cpp(Cp(patchi));

        htcBf[patchi] = 0.5*rhop*Cpp*magU*CfBf[patchi];
    }
}


// ************************************************************************* //

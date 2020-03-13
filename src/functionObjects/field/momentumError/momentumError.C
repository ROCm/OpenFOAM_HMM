/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "momentumError.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcLaplacian.H"
#include "turbulenceModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(momentumError, 0);
    addToRunTimeSelectionTable(functionObject, momentumError, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::functionObjects::momentumError::divDevRhoReff()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    {
        auto* turb = findObject<cmpTurbModel>
        (
            turbulenceModel::propertiesName
        );

        if (turb)
        {
            return tmp<volVectorField>::New
            (
                "divDevRhoReff",
              - fvc::div
                (
                    (turb->rho()*turb->nuEff())
                   *dev2(T(fvc::grad(turb->U()))),
                   "div(((rho*nuEff)*dev2(T(grad(U)))))"
                )
              - fvc::laplacian
                (
                    turb->rho()*turb->nuEff(),
                    turb->U(),
                    "laplacian(nuEff,U)"
                )
            );
        }
    }

    {
        const auto* turb = findObject<icoTurbModel>
        (
            turbulenceModel::propertiesName
        );

        if (turb)
        {
            return tmp<volVectorField>::New
            (
                "divDevReff",
              - fvc::div
                (
                    (turb->nuEff())*dev2(T(fvc::grad(turb->U()))),
                    "div((nuEff*dev2(T(grad(U)))))"
                )
              - fvc::laplacian
                (
                    turb->nuEff(), turb->U(), "laplacian(nuEff,U)"
                )
            );
        }
     }

     return volVectorField::null();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::momentumError::momentumError
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    pName_("p"),
    UName_("U"),
    phiName_("phi")
{
    read(dict);

    const surfaceScalarField& phi =
        lookupObject<surfaceScalarField>(phiName_);

    volVectorField* momentPtr
    (
        new volVectorField
        (
            IOobject
            (
                "momentError",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(phi.dimensions()*dimVelocity/dimVolume, Zero)
        )
    );

    mesh_.objectRegistry::store(momentPtr);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::momentumError::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Info<< type() << " " << name() << ":" << nl;

    // Optional field name entries
    if (dict.readIfPresent<word>("p", pName_))
    {
        Info<< "    p: " << pName_ << endl;
    }
    if (dict.readIfPresent<word>("U", UName_))
    {
        Info<< "    U: " << UName_ << endl;
    }

    if (dict.readIfPresent<word>("phi", phiName_))
    {
        Info<< "    phi: " << phiName_ << endl;
    }

    return true;
}


void Foam::functionObjects::momentumError::calcMomentError()
{

    volVectorField& momentErr =
        lookupObjectRef<volVectorField>("momentError");

    const volScalarField& p = lookupObject<volScalarField>(pName_);
    const volVectorField& U = lookupObject<volVectorField>(UName_);
    const surfaceScalarField& phi =
        lookupObject<surfaceScalarField>(phiName_);

    momentErr = divDevRhoReff() + fvc::div(phi, U) + fvc::grad(p);

}


bool Foam::functionObjects::momentumError::execute()
{
    calcMomentError();

    return true;
}


bool Foam::functionObjects::momentumError::write()
{
    const volVectorField& momentErr =
        lookupObjectRef<volVectorField>("momentError");

    momentErr.write();

    return true;
}


// ************************************************************************* //

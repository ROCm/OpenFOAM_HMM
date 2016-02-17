/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    applyBoundaryLayer

Description
    Apply a simplified boundary-layer model to the velocity and
    turbulence fields based on the 1/7th power-law.

    The uniform boundary-layer thickness is either provided via the -ybl option
    or calculated as the average of the distance to the wall scaled with
    the thickness coefficient supplied via the option -Cbl.  If both options
    are provided -ybl is used.

    Compressible modes is automatically selected based on the existence of the
    "thermophysicalProperties" dictionary required to construct the
    thermodynamics package.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallDist.H"
#include "processorFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// turbulence constants - file-scope
static const scalar Cmu(0.09);
static const scalar kappa(0.41);

void correctProcessorPatches(volScalarField& vf)
{
    if (!Pstream::parRun())
    {
        return;
    }

    // Not possible to use correctBoundaryConditions on fields as they may
    // use local info as opposed to the constraint values employed here,
    // but still need to update processor patches

    volScalarField::GeometricBoundaryField& bf = vf.boundaryField();

    forAll(bf, patchI)
    {
        if (isA<processorFvPatchField<scalar> >(bf[patchI]))
        {
            bf[patchI].initEvaluate();
        }
    }

    forAll(bf, patchI)
    {
        if (isA<processorFvPatchField<scalar> >(bf[patchI]))
        {
            bf[patchI].evaluate();
        }
    }
}


template<class TurbulenceModel>
Foam::tmp<Foam::volScalarField> calcK
(
    TurbulenceModel& turbulence,
    const volScalarField& mask,
    const volScalarField& nut,
    const volScalarField& y,
    const dimensionedScalar& ybl,
    const scalar Cmu,
    const scalar kappa
)
{
    // Turbulence k
    tmp<volScalarField> tk = turbulence->k();
    volScalarField& k = tk();
    scalar ck0 = pow025(Cmu)*kappa;
    k = (1 - mask)*k + mask*sqr(nut/(ck0*min(y, ybl)));
    k.rename("k");

    // Do not correct BC
    // - operation may use inconsistent fields wrt these local manipulations
    //k.correctBoundaryConditions();
    correctProcessorPatches(k);

    Info<< "Writing k\n" << endl;
    k.write();

    return tk;
}


template<class TurbulenceModel>
Foam::tmp<Foam::volScalarField> calcEpsilon
(
    TurbulenceModel& turbulence,
    const volScalarField& mask,
    const volScalarField& k,
    const volScalarField& y,
    const dimensionedScalar& ybl,
    const scalar Cmu,
    const scalar kappa
)
{
    // Turbulence epsilon
    tmp<volScalarField> tepsilon = turbulence->epsilon();
    volScalarField& epsilon = tepsilon();
    scalar ce0 = ::pow(Cmu, 0.75)/kappa;
    epsilon = (1 - mask)*epsilon + mask*ce0*k*sqrt(k)/min(y, ybl);
    epsilon.max(SMALL);
    epsilon.rename("epsilon");

    // Do not correct BC
    // - operation may use inconsistent fields wrt these local manipulations
    // epsilon.correctBoundaryConditions();
    correctProcessorPatches(epsilon);

    Info<< "Writing epsilon\n" << endl;
    epsilon.write();

    return tepsilon;
}


void calcOmega
(
    const fvMesh& mesh,
    const volScalarField& mask,
    const volScalarField& k,
    const volScalarField& epsilon
)
{
    // Turbulence omega
    IOobject omegaHeader
    (
        "omega",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (omegaHeader.headerOk())
    {
        volScalarField omega(omegaHeader, mesh);
        dimensionedScalar k0("SMALL", k.dimensions(), SMALL);

        omega = (1 - mask)*omega + mask*epsilon/(Cmu*k + k0);
        omega.max(SMALL);

        // Do not correct BC
        // - operation may use inconsistent fields wrt these local
        //   manipulations
        // omega.correctBoundaryConditions();
        correctProcessorPatches(omega);

        Info<< "Writing omega\n" << endl;
        omega.write();
    }
}


void setField
(
    const fvMesh& mesh,
    const word& fieldName,
    const volScalarField& value
)
{
    IOobject fldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (fldHeader.headerOk())
    {
        volScalarField fld(fldHeader, mesh);
        fld = value;

        // Do not correct BC
        // - operation may use inconsistent fields wrt these local
        //   manipulations
        // fld.correctBoundaryConditions();
        correctProcessorPatches(fld);

        Info<< "Writing " << fieldName << nl << endl;
        fld.write();
    }
}


void calcCompressible
(
    const fvMesh& mesh,
    const volScalarField& mask,
    const volVectorField& U,
    const volScalarField& y,
    const dimensionedScalar& ybl
)
{
    const Time& runTime = mesh.time();

    autoPtr<fluidThermo> pThermo(fluidThermo::New(mesh));
    fluidThermo& thermo = pThermo();
    volScalarField rho(thermo.rho());

    // Update/re-write phi
    #include "compressibleCreatePhi.H"
    phi.write();

    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    // Calculate nut - reference nut is calculated by the turbulence model
    // on its construction
    tmp<volScalarField> tnut = turbulence->nut();

    volScalarField& nut = tnut();
    volScalarField S(mag(dev(symm(fvc::grad(U)))));
    nut = (1 - mask)*nut + mask*sqr(kappa*min(y, ybl))*::sqrt(2)*S;

    // Do not correct BC - wall functions will 'undo' manipulation above
    // by using nut from turbulence model
    correctProcessorPatches(nut);
    nut.write();

    tmp<volScalarField> k =
        calcK(turbulence, mask, nut, y, ybl, Cmu, kappa);
    tmp<volScalarField> epsilon =
        calcEpsilon(turbulence, mask, k, y, ybl, Cmu, kappa);
    calcOmega(mesh, mask, k, epsilon);
    setField(mesh, "nuTilda", nut);
}


void calcIncompressible
(
    const fvMesh& mesh,
    const volScalarField& mask,
    const volVectorField& U,
    const volScalarField& y,
    const dimensionedScalar& ybl
)
{
    const Time& runTime = mesh.time();

    // Update/re-write phi
    #include "createPhi.H"
    phi.write();

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    tmp<volScalarField> tnut = turbulence->nut();

    // Calculate nut - reference nut is calculated by the turbulence model
    // on its construction
    volScalarField& nut = tnut();

    volScalarField S(mag(dev(symm(fvc::grad(U)))));
    nut = (1 - mask)*nut + mask*sqr(kappa*min(y, ybl))*::sqrt(2)*S;

    // Do not correct BC - wall functions will 'undo' manipulation above
    // by using nut from turbulence model
    correctProcessorPatches(nut);
    nut.write();

    tmp<volScalarField> k =
        calcK(turbulence, mask, nut, y, ybl, Cmu, kappa);
    tmp<volScalarField> epsilon =
        calcEpsilon(turbulence, mask, k, y, ybl, Cmu, kappa);
    calcOmega(mesh, mask, k, epsilon);
    setField(mesh, "nuTilda", nut);
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "apply a simplified boundary-layer model to the velocity and\n"
        "turbulence fields based on the 1/7th power-law."
    );

    #include "addRegionOption.H"

    argList::addOption
    (
        "ybl",
        "scalar",
        "specify the boundary-layer thickness"
    );
    argList::addOption
    (
        "Cbl",
        "scalar",
        "boundary-layer thickness as Cbl * mean distance to wall"
    );

    #include "setRootCase.H"

    if (!args.optionFound("ybl") && !args.optionFound("Cbl"))
    {
        FatalErrorInFunction
            << "Neither option 'ybl' or 'Cbl' have been provided to calculate "
            << "the boundary-layer thickness.\n"
            << "Please choose either 'ybl' OR 'Cbl'."
            << exit(FatalError);
    }
    else if (args.optionFound("ybl") && args.optionFound("Cbl"))
    {
        FatalErrorInFunction
            << "Both 'ybl' and 'Cbl' have been provided to calculate "
            << "the boundary-layer thickness.\n"
            << "Please choose either 'ybl' OR 'Cbl'."
            << exit(FatalError);
    }

    #include "createTime.H"
    #include "createNamedMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Modify velocity by applying a 1/7th power law boundary-layer
    // u/U0 = (y/ybl)^(1/7)
    // assumes U0 is the same as the current cell velocity

    Info<< "Setting boundary layer velocity" << nl << endl;
    scalar yblv = ybl.value();
    forAll(U, cellI)
    {
        if (y[cellI] <= yblv)
        {
            mask[cellI] = 1;
            U[cellI] *= ::pow(y[cellI]/yblv, (1.0/7.0));
        }
    }
    mask.correctBoundaryConditions();

    Info<< "Writing U\n" << endl;
    U.write();


    if
    (
        IOobject
        (
            basicThermo::dictName,
            runTime.constant(),
            mesh
        ).headerOk()
    )
    {
        calcCompressible(mesh, mask, U, y, ybl);
    }
    else
    {
        calcIncompressible(mesh, mask, U, y, ybl);
    }


    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

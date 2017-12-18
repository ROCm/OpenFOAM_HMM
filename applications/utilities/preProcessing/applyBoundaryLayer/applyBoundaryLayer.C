/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    applyBoundaryLayer

Group
    grpPreProcessingUtilities

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
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Turbulence constants - file-scope
static const scalar Cmu(0.09);
static const scalar kappa(0.41);


template<class Type>
void correctProcessorPatches
(
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    // Not possible to use correctBoundaryConditions on fields as they may
    // use local info as opposed to the constraint values employed here,
    // but still need to update processor patches
    typename volFieldType::Boundary& bf = vf.boundaryFieldRef();

    forAll(bf, patchi)
    {
        if (isA<processorFvPatchField<Type>>(bf[patchi]))
        {
            bf[patchi].initEvaluate();
        }
    }

    forAll(bf, patchi)
    {
        if (isA<processorFvPatchField<Type>>(bf[patchi]))
        {
            bf[patchi].evaluate();
        }
    }
}


void blendField
(
    const word& fieldName,
    const fvMesh& mesh,
    const scalarField& mask,
    const scalarField& boundaryLayerField
)
{
    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (fieldHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField fld(fieldHeader, mesh);
        scalarField& pf = fld.primitiveFieldRef();
        pf = (1 - mask)*pf + mask*boundaryLayerField;
        fld.max(SMALL);

        // Do not correct BC
        // - operation may use inconsistent fields wrt these local
        //   manipulations
        //fld.correctBoundaryConditions();
        correctProcessorPatches<scalar>(fld);

        Info<< "Writing " << fieldName << nl << endl;
        fld.write();
    }
}


void calcOmegaField
(
    const fvMesh& mesh,
    const scalarField& mask,
    const scalarField& kBL,
    const scalarField& epsilonBL
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

    if (omegaHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField omega(omegaHeader, mesh);
        scalarField& pf = omega.primitiveFieldRef();

        pf = (1 - mask)*pf + mask*epsilonBL/(Cmu*kBL + SMALL);
        omega.max(SMALL);

        // Do not correct BC
        // - operation may use inconsistent fields wrt these local
        //   manipulations
        // omega.correctBoundaryConditions();
        correctProcessorPatches<scalar>(omega);

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

    if (fldHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField fld(fldHeader, mesh);
        fld = value;

        // Do not correct BC
        // - operation may use inconsistent fields wrt these local
        //   manipulations
        // fld.correctBoundaryConditions();
        correctProcessorPatches<scalar>(fld);

        Info<< "Writing " << fieldName << nl << endl;
        fld.write();
    }
}


tmp<volScalarField> calcNut
(
    const fvMesh& mesh,
    const volVectorField& U
)
{
    const Time& runTime = mesh.time();

    if
    (
        IOobject
        (
            basicThermo::dictName,
            runTime.constant(),
            mesh
        ).typeHeaderOk<IOdictionary>(true)
    )
    {
        // Compressible
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

        // Correct nut
        turbulence->validate();

        return tmp<volScalarField>(new volScalarField(turbulence->nut()));
    }
    else
    {
        // Incompressible

        // Update/re-write phi
        #include "createPhi.H"
        phi.write();

        singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::turbulenceModel> turbulence
        (
            incompressible::turbulenceModel::New(U, phi, laminarTransport)
        );

        // Correct nut
        turbulence->validate();

        return tmp<volScalarField>(new volScalarField(turbulence->nut()));
    }
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
    argList::addBoolOption
    (
        "write-nut",
        "Write the turbulence viscosity field"
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

    bool writeNut = args.optionFound("write-nut");

    #include "createTime.H"
    #include "createNamedMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Modify velocity by applying a 1/7th power law boundary-layer
    // u/U0 = (y/ybl)^(1/7)
    // assumes U0 is the same as the current cell velocity
    Info<< "Setting boundary layer velocity" << nl << endl;
    scalar yblv = ybl.value();
    forAll(U, celli)
    {
        if (y[celli] <= yblv)
        {
            mask[celli] = 1;
            U[celli] *= ::pow(y[celli]/yblv, (1.0/7.0));
        }
    }
    mask.correctBoundaryConditions();
    correctProcessorPatches<vector>(U);

    // Retrieve nut from turbulence model
    volScalarField nut(calcNut(mesh, U));

    // Blend nut using boundary layer profile
    volScalarField S("S", mag(dev(symm(fvc::grad(U)))));
    nut = (1 - mask)*nut + mask*sqr(kappa*min(y, ybl))*::sqrt(2)*S;

    // Do not correct BC - wall functions will 'undo' manipulation above
    // by using nut from turbulence model
    correctProcessorPatches<scalar>(nut);
    if (writeNut)
    {
        Info<< "Writing nut\n" << endl;
        nut.write();
    }

    // Boundary layer turbulence kinetic energy
    scalar ck0 = pow025(Cmu)*kappa;
    scalarField kBL(sqr(nut/(ck0*min(y, ybl))));

    // Boundary layer turbulence dissipation
    scalar ce0 = ::pow(Cmu, 0.75)/kappa;
    scalarField epsilonBL(ce0*kBL*sqrt(kBL)/min(y, ybl));

    // Process fields if they are present
    blendField("k", mesh, mask, kBL);
    blendField("epsilon", mesh, mask, epsilonBL);
    calcOmegaField(mesh, mask, kBL, epsilonBL);
    if (writeNut) setField(mesh, "nuTilda", nut);

    // Write the updated U field
    Info<< "Writing U\n" << endl;
    U.write();

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

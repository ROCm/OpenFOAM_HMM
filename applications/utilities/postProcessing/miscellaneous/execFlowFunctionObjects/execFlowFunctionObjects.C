/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Global
    execFlowFunctionObjects

Description
    Execute the set of functionObjects specified in the selected dictionary
    (which defaults to system/controlDict) for the selected set of times.

    The flow (p-U) and optionally turbulence fields are available for the
    function objects to operate on allowing forces and other related properties
    to be calculated in addition to cutting planes etc.

\*---------------------------------------------------------------------------*/

#include "calc.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"

#include "incompressible/RASmodel/RASmodel.H"
#include "incompressible/LESmodel/LESmodel.H"

#include "basicThermo.H"
#include "compressible/RASmodel/RASmodel.H"
#include "compressible/LESmodel/LESmodel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    void execFlowFunctionObjects(const argList& args, const Time& runTime)
    {
        if (args.options().found("dict"))
        {
            fileName dictName(args.options()["dict"]);

            IOdictionary dict
            (
                IOobject
                (
                    dictName,
                    runTime.system(),
                    runTime,
                    IOobject::MUST_READ
                )
            );

            functionObjectList fol(runTime, dict);
            fol.start();
            fol.execute();
        }
        else
        {
            functionObjectList fol(runTime, runTime.controlDict());
            fol.start();
            fol.execute();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    Info<< "    Reading phi" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );

    Info<< "    Reading U" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );

    Info<< "    Reading p" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );

    if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        IOobject turbulencePropertiesHeader
        (
            "turbulenceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (turbulencePropertiesHeader.headerOk())
        {
            IOdictionary turbulenceProperties
            (
                turbulencePropertiesHeader
            );

            singlePhaseTransportModel laminarTransport(U, phi);

            if (turbulenceProperties.found("RASmodel"))
            {
                autoPtr<incompressible::RASmodel> RASmodel
                (
                    incompressible::RASmodel::New
                    (
                        U,
                        phi,
                        laminarTransport
                    )
                );

                execFlowFunctionObjects(args, runTime);
            }
            else if (turbulenceProperties.found("LESmodel"))
            {
                autoPtr<incompressible::LESmodel> sgsModel
                (
                    incompressible::LESmodel::New(U, phi, laminarTransport)
                );

                execFlowFunctionObjects(args, runTime);
            }
            else
            {
                FatalErrorIn(args.executable())
                    << "Cannot find turbulence model type in "
                    << "RASmodel dictionary"
                    << nl << exit(FatalError);
            }
        }
        else
        {
            IOdictionary transportProperties
            (
                IOobject
                (
                    "transportProperties",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            dimensionedScalar nu
            (
                transportProperties.lookup("nu")
            );

            execFlowFunctionObjects(args, runTime);
        }
    }
    else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
    {
        autoPtr<basicThermo> thermo
        (
            basicThermo::New(mesh)
        );

        volScalarField rho
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh
            ),
            thermo->rho()
        );

        IOobject turbulencePropertiesHeader
        (
            "turbulenceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (turbulencePropertiesHeader.headerOk())
        {
            IOdictionary turbulenceProperties
            (
                turbulencePropertiesHeader
            );

            if (turbulenceProperties.found("RASmodel"))
            {
                autoPtr<compressible::RASmodel> RASmodel
                (
                    compressible::RASmodel::New
                    (
                        rho,
                        U,
                        phi,
                        thermo()
                    )
                );

                execFlowFunctionObjects(args, runTime);
            }
            else if (turbulenceProperties.found("LESmodel"))
            {
                autoPtr<compressible::LESmodel> sgsModel
                (
                    compressible::LESmodel::New(rho, U, phi, thermo())
                );

                execFlowFunctionObjects(args, runTime);
            }
            else
            {
                FatalErrorIn(args.executable())
                    << "Cannot find turbulence model type in "
                    << "RASmodel dictionary"
                    << nl << exit(FatalError);
            }
        }
        else
        {
            IOdictionary transportProperties
            (
                IOobject
                (
                    "transportProperties",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            dimensionedScalar mu
            (
                transportProperties.lookup("mu")
            );

            execFlowFunctionObjects(args, runTime);
        }
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Incorrect dimensions of phi: " << phi.dimensions()
            << nl << exit(FatalError);
    }
}


// ************************************************************************* //

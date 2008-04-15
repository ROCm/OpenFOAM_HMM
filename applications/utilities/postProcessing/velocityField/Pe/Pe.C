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

Application
    Pe

Description
    Calculates and writes the Pe number as a surfaceScalarField obtained from
    field phi for each time

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"

#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "incompressible/LESmodel/LESmodel.H"

#include "basicThermo.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "compressible/LESmodel/LESmodel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject phiHeader
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check phi exists
        if (phiHeader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading phi" << endl;
            surfaceScalarField phi(phiHeader, mesh);

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


            Info<< "    Calculating Pe" << endl;

            if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
            {
                IOobject turbulencePropertiesHeader
                (
                    "turbulenceProperties",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (turbulencePropertiesHeader.headerOk())
                {
                    IOdictionary turbulenceProperties
                    (
                        turbulencePropertiesHeader
                    );

                    singlePhaseTransportModel laminarTransport(U, phi);

                    if (turbulenceProperties.found("turbulenceModel"))
                    {
                        autoPtr<turbulenceModel> turbulenceModel
                        (
                            turbulenceModel::New(U, phi, laminarTransport)
                        );

                        surfaceScalarField Pe
                        (
                            IOobject
                            (
                                "Pe",
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ
                            ),
                            mag(phi)
                            /(
                                mesh.magSf()
                               *mesh.surfaceInterpolation::deltaCoeffs()
                               *fvc::interpolate(turbulenceModel->nuEff())
                            )
                        );

                        Pe.write();
                    }
                    else if (turbulenceProperties.found("LESmodel"))
                    {
                        autoPtr<LESmodel> sgsModel
                        (
                            LESmodel::New(U, phi, laminarTransport)
                        );

                        surfaceScalarField Pe
                        (
                            IOobject
                            (
                                "Pe",
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ
                            ),
                            mag(phi)
                            /(
                                mesh.magSf()
                               *mesh.surfaceInterpolation::deltaCoeffs()
                               *fvc::interpolate(sgsModel->nuEff())
                            )
                        );

                        Info << "    Max Pe = " << max(Pe).value() << endl;

                        /*
                        label count = 0;
                        forAll(Pe, i)
                        {
                            if (Pe[i] > 200) count++;
                        }

                        Info<< "Fraction > 200 = "
                            << scalar(count)/Pe.size() << endl;
                        */

                        Pe.write();
                    }
                    else
                    {
                        FatalErrorIn(args.executable())
                            << "Cannot find turbulence model type in"
                               "turbulenceModel dictionary"
                            << exit(FatalError);
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

                    surfaceScalarField Pe
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mesh.surfaceInterpolation::deltaCoeffs()
                       *(mag(phi)/mesh.magSf())*(runTime.deltaT()/nu)
                    );

                    Info << "    Max Pe = " << max(Pe).value() << endl;
                    
                    Pe.write();
                }
            }
            else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
            {
                IOobject turbulencePropertiesHeader
                (
                    "turbulenceProperties",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (turbulencePropertiesHeader.headerOk())
                {
                    IOdictionary turbulenceProperties
                    (
                        turbulencePropertiesHeader
                    );

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

                    if (turbulenceProperties.found("turbulenceModel"))
                    {
                        autoPtr<compressible::turbulenceModel> turbulenceModel
                        (
                            compressible::turbulenceModel::New
                            (
                                rho,
                                U,
                                phi,
                                thermo()
                            )
                        );

                        surfaceScalarField Pe
                        (
                            IOobject
                            (
                                "Pe",
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ
                            ),
                            mag(phi)
                            /(
                                mesh.magSf()
                               *mesh.surfaceInterpolation::deltaCoeffs()
                               *fvc::interpolate(turbulenceModel->muEff())
                            )
                        );

                        Pe.write();
                    }
                    else if (turbulenceProperties.found("LESmodel"))
                    {
                        autoPtr<compressible::LESmodel> sgsModel
                        (
                            compressible::LESmodel::New(rho, U, phi, thermo())
                        );

                        surfaceScalarField Pe
                        (
                            IOobject
                            (
                                "Pe",
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ
                            ),
                            mag(phi)
                            /(
                                mesh.magSf()
                               *mesh.surfaceInterpolation::deltaCoeffs()
                               *fvc::interpolate(sgsModel->muEff())
                            )
                        );

                        Info << "    Max Pe = " << max(Pe).value() << endl;

                        Pe.write();
                    }
                    else
                    {
                        FatalErrorIn(args.executable())
                            << "Cannot find turbulence model type in"
                               "turbulenceModel dictionary"
                            << exit(FatalError);
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
 
                    surfaceScalarField Pe
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mesh.surfaceInterpolation::deltaCoeffs()
                       *(mag(phi)/(mesh.magSf()))*(runTime.deltaT()/mu)
                    );

                    Info << "    Max Pe = " << max(Pe).value() << endl;
                
                    Pe.write();
                }
            }
            else
            {
                FatalErrorIn(args.executable())
                    << "Incorrect dimensions of phi: " << phi.dimensions()
                    << abort(FatalError);
            }
        }
        else
        {
            Info<< "    No phi" << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

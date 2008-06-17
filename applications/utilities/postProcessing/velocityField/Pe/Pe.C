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
    field phi.
    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RASmodel/RASmodel.H"
#include "incompressible/LESmodel/LESmodel.H"
#include "basicThermo.H"
#include "compressible/RASmodel/RASmodel.H"
#include "compressible/LESmodel/LESmodel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.options().found("noWrite");

    IOobject phiHeader
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (phiHeader.headerOk())
    {
        autoPtr<surfaceScalarField> PePtr;

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

                    PePtr.set
                    (
                        new surfaceScalarField
                        (
                            IOobject
                            (
                                "Pe",
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ
                            ),
                            mag(phi) /
                            (
                                mesh.magSf()
                              * mesh.surfaceInterpolation::deltaCoeffs()
                              * fvc::interpolate(RASmodel->nuEff())
                            )
                        )
                    );
                }
                else if (turbulenceProperties.found("LESmodel"))
                {
                    autoPtr<LESmodel> sgsModel
                    (
                        LESmodel::New(U, phi, laminarTransport)
                    );

                    PePtr.set
                    (
                        new surfaceScalarField
                        (
                            IOobject
                            (
                                "Pe",
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ
                            ),
                            mag(phi) /
                            (
                                mesh.magSf()
                              * mesh.surfaceInterpolation::deltaCoeffs()
                              * fvc::interpolate(sgsModel->nuEff())
                            )
                        )
                    );
                }
                else
                {
                    FatalErrorIn(args.executable())
                        << "Cannot find turbulence model type in "
                        "RASmodel dictionary"
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

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mesh.surfaceInterpolation::deltaCoeffs()
                      * (mag(phi)/mesh.magSf())*(runTime.deltaT()/nu)
                    )
                );
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

                    PePtr.set
                    (
                        new surfaceScalarField
                        (
                            IOobject
                            (
                                "Pe",
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ
                            ),
                            mag(phi) /
                            (
                                mesh.magSf()
                              * mesh.surfaceInterpolation::deltaCoeffs()
                              * fvc::interpolate(RASmodel->muEff())
                            )
                        )
                    );
                }
                else if (turbulenceProperties.found("LESmodel"))
                {
                    autoPtr<compressible::LESmodel> sgsModel
                    (
                        compressible::LESmodel::New(rho, U, phi, thermo())
                    );

                    PePtr.set
                    (
                        new surfaceScalarField
                        (
                            IOobject
                            (
                                "Pe",
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ
                            ),
                            mag(phi) /
                            (
                                mesh.magSf()
                              * mesh.surfaceInterpolation::deltaCoeffs()
                              * fvc::interpolate(sgsModel->muEff())
                            )
                        )
                    );
                }
                else
                {
                    FatalErrorIn(args.executable())
                        << "Cannot find turbulence model type in"
                        "RASmodel dictionary"
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

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mesh.surfaceInterpolation::deltaCoeffs()
                      * (mag(phi)/(mesh.magSf()))*(runTime.deltaT()/mu)
                    )
                );
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Incorrect dimensions of phi: " << phi.dimensions()
                    << abort(FatalError);
        }


        // can also check how many cells exceed a particular Pe limit
        /*
        {
            label count = 0;
            label PeLimit = 200;
            forAll(PePtr(), i)
            {
                if (PePtr()[i] > PeLimit)
                {
                    count++;
                }

            }

            Info<< "Fraction > " << PeLimit << " = "
                << scalar(count)/Pe.size() << endl;
        }
        */

        Info << "Pe max : " << max(PePtr()).value() << endl;

        if (writeResults)
        {
            PePtr().write();
        }
    }
    else
    {
        Info<< "    No phi" << endl;
    }
}

// ************************************************************************* //

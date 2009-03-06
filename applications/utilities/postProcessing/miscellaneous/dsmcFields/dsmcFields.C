/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    dsmcFields

Description
    Calculate intensive fields (U and T) from averaged extensive fields from a
    DSMC calculation.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.options().found("noWrite");

    IOobject rhoNMeanheader
    (
        "rhoNMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject rhoMMeanheader
    (
        "rhoMMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject momentumMeanheader
    (
        "momentumMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject linearKEMeanheader
    (
        "linearKEMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject internalEMeanheader
    (
        "internalEMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject iDofMeanheader
    (
        "iDofMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (!rhoNMeanheader.headerOk())
    {
        Info<< "    No rhoNMean" << endl;
    }
    else if (!rhoMMeanheader.headerOk())
    {
        Info<< "    No rhoMMean" << endl;
    }
    else if (!momentumMeanheader.headerOk())
    {
        Info<< "    No momentumMean" << endl;
    }
    else if (!linearKEMeanheader.headerOk())
    {
        Info<< "    No linearKEMean" << endl;
    }
    else if (!internalEMeanheader.headerOk())
    {
        Info<< "    No internalEMean" << endl;
    }
    else if (!iDofMeanheader.headerOk())
    {
        Info<< "    No iDofMean" << endl;
    }
    else
    {
        Info<< "Reading field rhoNMean" << endl;
        volScalarField rhoNMean(rhoNMeanheader, mesh);

        Info<< "Reading field rhoMMean" << endl;
        volScalarField rhoMMean(rhoMMeanheader, mesh);

        Info<< "Reading field momentumMean" << endl;
        volVectorField momentumMean(momentumMeanheader, mesh);

        Info<< "Reading field linearKEMean" << endl;
        volScalarField linearKEMean(linearKEMeanheader, mesh);

        Info<< "Reading field internalEMean" << endl;
        volScalarField internalEMean(internalEMeanheader, mesh);

        Info<< "Reading field iDofMean" << endl;
        volScalarField iDofMean(iDofMeanheader, mesh);

        // Check if there are any zero values in the density fields

        if (min(rhoNMean).value() > VSMALL)
        {
            // TODO Sort out boundary field values if required

            Info<< nl << "Calculating UMean field." << endl;
            volVectorField UMean
            (
                IOobject
                (
                    "UMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                momentumMean/rhoMMean
            );

            Info<< nl << "Calculating translationalTMean field." << endl;
            volScalarField translationalTMean
            (
                IOobject
                (
                    "translationalTMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                2.0/(3.0*dsmcCloud::kb*rhoNMean)
               *(linearKEMean - 0.5*rhoMMean*(UMean & UMean))
            );

            Info<< nl << "Calculating internalTMean field." << endl;
            volScalarField internalTMean
            (
                IOobject
                (
                    "internalTMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                2.0/(dsmcCloud::kb*iDofMean)*internalEMean
            );

            Info<< nl << "Calculating overallTMean field." << endl;
            volScalarField overallTMean
            (
                IOobject
                (
                    "overallTMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                2.0/(dsmcCloud::kb*(3.0*rhoNMean + iDofMean))
               *(linearKEMean - 0.5*rhoMMean*(UMean & UMean) + internalEMean)
            );

            Info<< nl << "mag(UMean) max/min : "
                << max(mag(UMean)).value() << " "
                << min(mag(UMean)).value() << endl;

            Info<< nl << "translationalTMean max/min : "
                << max(translationalTMean).value() << " "
                << min(translationalTMean).value() << endl;

            Info<< nl << "internalTMean max/min : "
                << max(internalTMean).value() << " "
                << min(internalTMean).value() << endl;

            Info<< nl << "overallTMean max/min : "
                << max(overallTMean).value() << " "
                << min(overallTMean).value() << endl;

            if (writeResults)
            {
                UMean.write();
                translationalTMean.write();
                internalTMean.write();
                overallTMean.write();
            }
        }
        else
        {
            Info<< "Small or negative value (" << min(rhoNMean)
                << ") found in rhoNMean field. "
                << "Not calculating fields to avoid division by zero "
                << "or invalid results."
                << endl;
        }
    }
}

// ************************************************************************* //

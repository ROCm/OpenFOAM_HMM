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
    flowType

Description
    Calculates and writes the flowType of velocity field U at each time

    The flow tye parameter is obtained from

                 |D| - |Omega|
        lambda = -------------
                 |D| + |Omega|

        -1 = rotational flow
         0 = simple shear flow
         1 = planar extensional flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

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

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check U exists
        if (Uheader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

            volTensorField gradU = fvc::grad(U);
            volScalarField magD = mag(symm(gradU));
            volScalarField magOmega = mag(skew(gradU));
            dimensionedScalar smallMagD("smallMagD", magD.dimensions(), SMALL);

            Info<< "    Calculating flowType" << endl;
            volScalarField flowType
            (
                IOobject
                (
                    "flowType",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                (magD - magOmega)/(magD + magOmega + smallMagD)
            );

            flowType.write();
        }
        else
        {
            Info<< "    No U" << endl;
        }
    }

    return(0);
}


// ************************************************************************* //

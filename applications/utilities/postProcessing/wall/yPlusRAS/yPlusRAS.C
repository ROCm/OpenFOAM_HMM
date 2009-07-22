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
    yPlusRAS

Description
    Calculates and reports yPlus for all wall patches, for the specified times.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "RASModel.H"
#include "wallFvPatch.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        fvMesh::readUpdateState state = mesh.readUpdate();

        // Wall distance
        if (timeI == 0 || state != fvMesh::UNCHANGED)
        {
            Info<< "Calculating wall distance\n" << endl;
            wallDist y(mesh, true);
            Info<< "Writing wall distance to field "
                << y.name() << nl << endl;
            y.write();
        }

        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );

        Info << "Reading field U\n" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

#       include "createPhi.H"

        singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::RASModel> rasModel
        (
            incompressible::RASModel::New(U, phi, laminarTransport)
        );

        const scalar Cmu =
            rasModel->coeffDict().lookupOrDefault<scalar>("Cmu", 0.09);

        const fvPatchList& patches = mesh.boundary();

        forAll(patches, patchI)
        {
            const fvPatch& currPatch = patches[patchI];

            if (typeid(currPatch) == typeid(wallFvPatch))
            {
                yPlus.boundaryField()[patchI] = rasModel->yPlus(patchI, Cmu);
                const scalarField& Yp = yPlus.boundaryField()[patchI];

                Info<< "Patch " << patchI
                    << " named " << currPatch.name()
                    << " y+ : min: " << min(Yp) << " max: " << max(Yp)
                    << " average: " << average(Yp) << nl << endl;
            }
        }

        Info<< "Writing yPlus to field "
            << yPlus.name() << nl << endl;

        yPlus.write();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //

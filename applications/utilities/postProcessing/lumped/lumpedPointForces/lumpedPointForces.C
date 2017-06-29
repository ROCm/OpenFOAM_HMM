/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

Application
    lumpedPointForces

Description
    Extract force/moment information from existing calculations based
    on the segmentation of the pressure integration zones.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "volFields.H"
#include "IOobjectList.H"

#include "lumpedPointTools.H"
#include "lumpedPointIOMovement.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoFieldType>
autoPtr<GeoFieldType> loadField
(
    const volMesh::Mesh& mesh,
    const IOobject* io
)
{
    if (io && io->headerClassName() == GeoFieldType::typeName)
    {
        Info<< "Reading " << GeoFieldType::typeName
            << ' ' << io->name() << endl;

        return autoPtr<GeoFieldType>
        (
            new GeoFieldType
            (
                IOobject
                (
                    io->name(),
                    io->instance(),
                    io->local(),
                    io->db(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE,
                    io->registerObject()
                ),
                mesh
            )
        );
    }

    return autoPtr<GeoFieldType>();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Extract force/moment information from existing calculations based "
        "on the lumpedPoints pressure zones."
    );

    argList::addBoolOption
    (
        "vtk",
        "create visualization files of the forces"
    );

    timeSelector::addOptions(true, false);
    argList::noFunctionObjects();  // Never use function objects

    #include "addRegionOption.H"
    #include "setRootCase.H"

    const bool withVTK = args.optionFound("vtk");

    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    autoPtr<lumpedPointIOMovement> movement = lumpedPointIOMovement::New
    (
        runTime
    );

    if (!movement.valid())
    {
        Info<< "no valid movement given" << endl;
        return 1;
    }

    const labelList patchLst = lumpedPointTools::lumpedPointPatchList(mesh);
    if (patchLst.empty())
    {
        Info<< "no patch list found" << endl;
        return 2;
    }

    movement().setMapping(mesh, patchLst, lumpedPointTools::points0Field(mesh));

    List<vector> forces, moments;
    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        if (mesh.readUpdate())
        {
            Info<< "    Read new mesh" << nl;
        }

        // Search for list of objects for this time
        IOobjectList objects(mesh, runTime.timeName());

        // Pressure field
        autoPtr<volScalarField> pField
            = loadField<volScalarField>(mesh, objects.lookup("p"));

        // The forces per zone
        if (movement().forcesAndMoments(mesh, forces, moments))
        {
            Info<<"forces per zone: " << forces << endl;
            Info<<"moments per zone: " << moments << endl;

            if (withVTK && Pstream::master())
            {
                const word outputName =
                    Foam::name("forces_%06d.vtp", runTime.timeIndex());

                Info<<"    " << outputName << endl;

                movement().writeForcesAndMomentsVTP
                (
                    outputName,
                    forces,
                    moments
                );
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //

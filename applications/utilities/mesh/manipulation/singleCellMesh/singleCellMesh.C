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
    singleCellMesh

Description
    Removes all but one cells of the mesh. Used to generate mesh and fields
    that can be used for boundary-only data.
    Might easily result in illegal mesh though so only look at boundaries
    in paraview.

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Time.H"
#include "ReadFields.H"
#include "singleCellFvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void interpolateFields
(
    const singleCellFvMesh& scMesh,
    const PtrList<GeoField>& flds
)
{
    forAll(flds, i)
    {
        tmp<GeoField> scFld = scMesh.interpolate(flds[i]);
        GeoField* scFldPtr = scFld.ptr();
        scFldPtr->writeOpt() = IOobject::AUTO_WRITE;
        scFldPtr->store();
    }
}


// Main program:

int main(int argc, char *argv[])
{
#   include "addOverwriteOption.H"
#   include "addTimeOptions.H"

#   include "setRootCase.H"
#   include "createTime.H"
    // Get times list
    instantList Times = runTime.times();
#   include "checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);
#   include "createMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const bool overwrite = args.optionFound("overwrite");


    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.
    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);


    if (!overwrite)
    {
        runTime++;
    }

    // Create the mesh
    singleCellFvMesh scMesh
    (
        IOobject
        (
            mesh.name(),
            mesh.polyMesh::instance(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    // Map and store the fields on the scMesh.
    interpolateFields(scMesh, vsFlds);
    interpolateFields(scMesh, vvFlds);
    interpolateFields(scMesh, vstFlds);
    interpolateFields(scMesh, vsymtFlds);
    interpolateFields(scMesh, vtFlds);


    // Write
    Info<< "Writing mesh to time " << runTime.timeName() << endl;
    scMesh.write();


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

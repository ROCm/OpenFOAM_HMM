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
    applyWallFunctionBounaryConditions

Description
    Updates OpenFOAM incompressible RAS cases to use the new wall function
    framework

    NOTE: For incompressible RAS calculations ONLY

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"

#include "wallPolyPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main program:

void createNut(const fvMesh& mesh)
{
    IOobject nutHeader
    (
        "nut",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!nutHeader.headerOk())
    {
        Info<< "Creating field nut" << nl << endl;

        volScalarField nut
        (
            IOobject
            (
                "nut",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
        );

        nut.write();
    }
}


void replaceBoundaryType
(
    const fvMesh& mesh,
    const word& fieldName,
    const word& boundaryType,
    const string& boundaryValue
)
{
    IOobject header
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!header.headerOk())
    {
        return;
    }

    Info<< "Updating boundary types for field " << header.name() << endl;

    const word oldTypeName = IOdictionary::typeName;
    const_cast<word&>(IOdictionary::typeName) = word::null;

    IOdictionary dict(header);

    const_cast<word&>(IOdictionary::typeName) = oldTypeName;
    const_cast<word&>(dict.type()) = dict.headerClassName();

    // Make a backup of the old field
    word backupName(dict.name() + ".old");
    Info<< "    copying original " << dict.name() << " to "
        << backupName << endl;
    IOdictionary dictOld = dict;
    dictOld.rename(backupName);
    dictOld.regIOobject::write();

    // Loop through boundary patches and update
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    dictionary& boundaryDict = dict.subDict("boundaryField");
    forAll(bMesh, patchI)
    {
        if (isType<wallPolyPatch>(bMesh[patchI]))
        {
            word patchName = bMesh[patchI].name();
            dictionary& oldPatch = boundaryDict.subDict(patchName);

            dictionary newPatch(dictionary::null);
            newPatch.add("type", boundaryType);
            newPatch.add("value", ("uniform " + boundaryValue).c_str());

            oldPatch = newPatch;
        }
    }

    Info<< "    writing updated " << dict.name() << nl << endl;
    dict.regIOobject::write();
}


int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    Info<< "Updating turbulence fields to operate using new run time "
        << "selectable" << nl << "wall functions" << nl << nl
        << ">>>>NOTE: only applicable to incompressible RAS models"
        << nl << endl;

    createNut(mesh);

    replaceBoundaryType(mesh, "nut", "nutWallFunction", "0");
    replaceBoundaryType(mesh, "epsilon", "epsilonWallFunction", "0");
    replaceBoundaryType(mesh, "omega", "omegaWallFunction", "0");
    replaceBoundaryType(mesh, "k", "kQRWallFunction", "0");
    replaceBoundaryType(mesh, "q", "kQRWallFunction", "0");
    replaceBoundaryType(mesh, "R", "kQRWallFunction", "(0 0 0 0 0 0)");

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

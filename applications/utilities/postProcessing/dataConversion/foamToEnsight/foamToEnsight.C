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

Description
    Translates FOAM data to EnSight format

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"
#include "tensorIOField.H"
#include "IOobjectList.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "scalarIOField.H"

#include "ensightMesh.H"
#include "ensightField.H"

#include "ensightParticlePositions.H"
#include "ensightSprayField.H"

#include "fvc.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("patches", "patch list");
    argList::validOptions.insert("binary", "" );
#   include "addTimeOptions.H"

    const label nTypes = 2;
    const word fieldTypes[] =
    {
        volScalarField::typeName,
        volVectorField::typeName
    };

    const label nSprayFieldTypes = 2;
    const word sprayFieldTypes[] =
    {
        scalarIOField::typeName,
        vectorIOField::typeName
    };

#   include "setRootCase.H"
#   include "createTime.H"

    // get the available time-steps
    instantList Times = runTime.times();

#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"

    const word postProcDir = "EnSight";
    const word prepend = args.globalCaseName() + '.';
    const word sprayName = "lagrangian";

    fileName postProcPath = args.rootPath()/args.globalCaseName()/postProcDir;

    if (Pstream::master())
    {
        if (dir(postProcPath))
        {
            rmDir(postProcPath);
        }

        mkDir(postProcPath);
    }

    OFstream *ensightCaseFilePtr = NULL;

    // Check options
    bool binary = false;
    if (args.options().found("binary"))
    {
        binary = true;
    }

    if (Pstream::master())
    {
        // Open the Case file
        fileName ensightCaseFileName = prepend + "case";

        if (!binary)
        {
            ensightCaseFilePtr = new OFstream
            (
                postProcPath/ensightCaseFileName,
                runTime.writeFormat(),
                runTime.writeVersion(),
                runTime.writeCompression()
            );
        }
        else
        {
            ensightCaseFilePtr = new OFstream
            (
                postProcPath/ensightCaseFileName,
                runTime.writeFormat(),
                runTime.writeVersion(),
                IOstream::UNCOMPRESSED
            );
        }

        Info<< nl << "Case file is " << ensightCaseFileName << endl;
    }

    OFstream& ensightCaseFile = *ensightCaseFilePtr;

    // Construct the EnSight mesh
    ensightMesh eMesh(mesh, args, binary);

    // Set Time to the last time before looking for the spray objects
    runTime.setTime(Times[Times.size()-1], Times.size()-1);

    IOobjectList objects(mesh, runTime.timeName());
    IOobjectList sprayObjects(mesh, runTime.timeName(), "lagrangian");

    bool lagrangianExist = false;

    if (!eMesh.patchNames.size())
    {
        IOobject lagrangianHeader
        (
            "positions",
            runTime.timeName(),
            "lagrangian",
            mesh,
            IOobject::NO_READ
        );

        if (lagrangianHeader.headerOk())
        {
            lagrangianExist = true;
        }
    }

#   include "ensightCaseHeader.H"

#   include "checkMeshMoving.H"

    word geomCaseFileName = prepend + "000";
    if (Pstream::master())
    {
        // test pre check variable if there is a moving mesh
        if (meshMoving == true) geomCaseFileName = prepend + "***";
        ensightCaseFile
            << "GEOMETRY" << nl
            << "model:        1     "
            << (geomCaseFileName + ".mesh").c_str() << nl;
    }

    label nTimeSteps = 0;
    for (label n=startTime; n<endTime; n++)
    {
        nTimeSteps++;
        runTime.setTime(Times[n], n);
        label timeIndex = n - startTime;

        word timeName = itoa(timeIndex);
        word timeFile = prepend + timeName;

        Info << "Translating time = " << runTime.timeName() << nl;

#       include "moveMesh.H"

        if (timeIndex == 0 || mesh.moving())
        {
            eMesh.write
            (
                postProcPath,
                prepend,
                timeIndex,
                ensightCaseFile
            );
        }

        if (Pstream::master() && timeIndex == 0)
        {
            if (lagrangianExist)
            {
                ensightCaseFile
                    <<  (
                            "measured:     1     "
                          + prepend
                          + "***."
                          + sprayName
                        ).c_str()
                    << nl;
             }
             ensightCaseFile << nl << "VARIABLE" << nl;
        }

        for (label i=0; i<nTypes; i++)
        {
            wordList fieldNames = objects.names(fieldTypes[i]);

            for (label j=0; j<fieldNames.size(); j++)
            {
                word fieldName = fieldNames[j];

                bool variableGood = true;

#               include "checkData.H"

                if (variableGood)
                {
                    IOobject fieldObject
                    (
                        fieldName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    if (fieldTypes[i] == volScalarField::typeName)
                    {
                        ensightField<scalar>
                        (
                            fieldObject,
                            eMesh,
                            postProcPath,
                            prepend,
                            timeIndex,
                            binary,
                            ensightCaseFile
                        );
                    }
                    else if (fieldTypes[i] == volVectorField::typeName)
                    {
                        ensightField<vector>
                        (
                            fieldObject,
                            eMesh,
                            postProcPath,
                            prepend,
                            timeIndex,
                            binary,
                            ensightCaseFile
                        );
                    }
                    else if (fieldTypes[i] == volSphericalTensorField::typeName)
                    {
                        ensightField<sphericalTensor>
                        (
                            fieldObject,
                            eMesh,
                            postProcPath,
                            prepend,
                            timeIndex,
                            binary,
                            ensightCaseFile
                        );
                    }
                    else if (fieldTypes[i] == volSymmTensorField::typeName)
                    {
                        ensightField<symmTensor>
                        (
                            fieldObject,
                            eMesh,
                            postProcPath,
                            prepend,
                            timeIndex,
                            binary,
                            ensightCaseFile
                        );
                    }
                    else if (fieldTypes[i] == volTensorField::typeName)
                    {
                        ensightField<tensor>
                        (
                            fieldObject,
                            eMesh,
                            postProcPath,
                            prepend,
                            timeIndex,
                            binary,
                            ensightCaseFile
                        );
                    }
                }
            }
        }

        if (lagrangianExist)
        {
            ensightParticlePositions
            (
                mesh,
                postProcPath,
                timeFile,
                sprayName
            );

            for (label i=0; i<nSprayFieldTypes; i++)
            {
                wordList fieldNames = sprayObjects.names(sprayFieldTypes[i]);

                for (label j=0; j<fieldNames.size(); j++)
                {
                    word fieldName = fieldNames[j];

                    IOobject fieldObject
                    (
                        fieldName,
                        mesh.time().timeName(),
                        "lagrangian",
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    if (sprayFieldTypes[i] == scalarIOField::typeName)
                    {
                        ensightSprayField<scalar>
                        (
                            fieldObject,
                            postProcPath,
                            prepend,
                            timeIndex,
                            sprayName,
                            ensightCaseFile
                        );
                    }
                    else if (sprayFieldTypes[i] == vectorIOField::typeName)
                    {
                        ensightSprayField<vector>
                        (
                            fieldObject,
                            postProcPath,
                            prepend,
                            timeIndex,
                            sprayName,
                            ensightCaseFile
                        );
                    }
                    else if (sprayFieldTypes[i] == tensorIOField::typeName)
                    {
                        ensightSprayField<tensor>
                        (
                            fieldObject,
                            postProcPath,
                            prepend,
                            timeIndex,
                            sprayName,
                            ensightCaseFile
                        );
                    }
                }
            }
        }
    }

#   include "ensightCaseTail.H"

    if (Pstream::master())
    {
        delete ensightCaseFilePtr;
    }

    return 0;
}


// ************************************************************************* //

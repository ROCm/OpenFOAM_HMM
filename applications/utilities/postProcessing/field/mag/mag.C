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
    mag

Description
    Calculates and writes the magnitude of a field for each time

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "writeMagField.C"

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("fieldName1 .. fieldNameN"); // abuse for usage

    // setRootCase, but skip args check
    argList args(argc, argv, false);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

    const stringList& params = args.additionalArgs();
    if (!params.size())
    {
        Info<< nl << "must specify one or more fields" << nl;
        args.printUsage();
        FatalError.exit();
    }

#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        forAll(params, paramI)
        {
            const word fieldName(params[paramI]);

            IOobject fieldHeader
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );

            // Check field exists
            if (fieldHeader.headerOk())
            {
                bool processed = false;

                writeMagField<scalar>(fieldHeader, mesh, processed);
                writeMagField<vector>(fieldHeader, mesh, processed);
                writeMagField<sphericalTensor>(fieldHeader, mesh, processed);
                writeMagField<symmTensor>(fieldHeader, mesh, processed);
                writeMagField<tensor>(fieldHeader, mesh, processed);

                if (!processed)
                {
                    FatalError
                        << "Unable to process " << fieldName << nl
                        << "No call to mag for fields of type "
                        << fieldHeader.headerClassName() << nl << nl
                        << exit(FatalError);
                }
            }
            else
            {
                Info<< "    No " << fieldName << endl;
            }
        }
    }

    return 0;
}

// ************************************************************************* //

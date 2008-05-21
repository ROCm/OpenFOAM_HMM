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
    magGrad

Description
    Calculates and writes the scalar magnitude of a scalar or vector field
    at each time

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("field1 ... fieldN");   // abuse for usage

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
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        forAll(params, paramI)
        {
            const word fieldName(params[paramI]);

            IOobject header
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );

            // Check field exists
            if (header.headerOk())
            {
                if (header.headerClassName() == volScalarField::typeName)
                {
                    Info<< "    Reading " << fieldName << " ...";
                    volScalarField field(header, mesh);

                    volScalarField magGrad
                    (
                        IOobject
                        (
                            "magGrad" + fieldName,
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(fvc::grad(field))
                    );
                    Info<< "Writing " << magGrad.name() << endl;
                    magGrad.write();
                }
                else if (header.headerClassName() == volVectorField::typeName)
                {
                    Info<< "    Reading " << fieldName << " ...";
                    volVectorField field(header, mesh);

                    volScalarField magGrad
                    (
                        IOobject
                        (
                            "magGrad" + fieldName,
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(fvc::grad(field))
                    );
                    Info<< "Writing " << magGrad.name() << endl;
                    magGrad.write();
                }
                else
                {
                    Info<< "    Skipping " << fieldName << " : type "
                        << header.headerClassName() << endl;
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

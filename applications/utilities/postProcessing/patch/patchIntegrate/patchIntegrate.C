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
    patchIntegrate

Description
    Calculates the integral of the specified field over the specified patch.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("fieldName");
    argList::validArgs.append("patchName");

#   include "addTimeOptions.H"
#   include "setRootCase.H"

    word fieldName(args.additionalArgs()[0]);
    word patchName(args.additionalArgs()[1]);

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
            mesh.readUpdate();

            label patchi = mesh.boundaryMesh().findPatchID(patchName);
            if (patchi < 0)
            {
                FatalError
                    << "Unable to find patch " << patchName << nl
                    << exit(FatalError);
            }

            if (fieldHeader.headerClassName() == "volScalarField")
            {
                Info<< "    Reading volScalarField " << fieldName << endl;
                volScalarField field(fieldHeader, mesh);

                vector sumField = sum
                    (
                        mesh.Sf().boundaryField()[patchi]
                       *field.boundaryField()[patchi]
                    );

                Info<< "    Integral of " << fieldName << " over patch "
                    << patchName << '[' << patchi << ']' << " = "
                    << sumField << nl;
            }
            else if (fieldHeader.headerClassName() == "surfaceScalarField")
            {
                Info<< "    Reading surfaceScalarField " << fieldName << endl;

                surfaceScalarField field(fieldHeader, mesh);
                scalar sumField = sum(field.boundaryField()[patchi]);

                Info<< "    Integral of " << fieldName << " over patch "
                    << patchName << '[' << patchi << ']' << " = "
                    << sumField << nl;
            }
            else
            {
                FatalError
                    << "Only possible to integrate volScalarFields "
                    << "and surfaceScalarFields" << nl << exit(FatalError);
            }
        }
        else
        {
            Info<< "    No field " << fieldName << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

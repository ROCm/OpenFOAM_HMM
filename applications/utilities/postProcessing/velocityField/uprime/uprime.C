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
    uprime

Description
    Calculates and writes the scalar field of uprime (sqrt(2/3 k)) at 
    each time

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

        IOobject kheader
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check k exists
        if (kheader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading k" << endl;
            volScalarField k(kheader, mesh);

            Info<< "    Calculating uprime" << endl;
            volScalarField uprime
            (
                IOobject
                (
                    "uprime",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                sqrt((2.0/3.0)*k)
            );
            uprime.write();
        }
        else
        {
            Info<< "    No k" << endl;
        }

        Info<< endl;
    }

    return(0);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    star3ToFoam

Description
    Converts a Star-CD (v3) pro-STAR mesh into OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "starMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("STAR mesh file prefix");
    argList::addOption("scale", "scale factor");

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    scalar scaleFactor = args.optionLookupOrDefault("scale", 1.0);

#   include "createTime.H"

    starMesh makeMesh(args[1], runTime, scaleFactor);

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info<< "Writing mesh" << endl;
    makeMesh.writeMesh();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

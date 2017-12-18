/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    fireToFoam

Group
    grpMeshConversionUtilities

Description
    Converts an AVL/FIRE polyhedral mesh to OPENFOAM

Usage
    \b fireToFoam [OPTION] firePolyMesh

    Options:

      - \par -ascii
        Write in ASCII format instead of binary

      - \par -check
        Perform edge checking

      - \par -scale \<factor\>
        Specify an alternative geometry scaling factor.
        The default is \b 1 (no scaling).

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "FIREMeshReader.H"
#include "checkFireEdges.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert AVL/FIRE polyhedral mesh to OPENFOAM format"
    );

    argList::noParallel();
    argList::addArgument("firePolyMesh");
    argList::addBoolOption
    (
        "ascii",
        "write in ASCII format instead of binary"
    );
    argList::addBoolOption
    (
        "check",
        "perform edge checking as well"
    );
    argList::addOption
    (
        "scale",
        "scale",
        "geometry scaling factor - default is 1 (no scaling)"
    );


    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());


    // Binary output, unless otherwise specified
    const IOstream::streamFormat format =
    (
        args.optionFound("ascii")
      ? IOstream::ASCII
      : IOstream::BINARY
    );

    // increase the precision of the points data
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));


    fileFormats::FIREMeshReader reader
    (
        args[1],
        // Default no scaling
        args.optionLookupOrDefault("scale", 1.0)
    );


    autoPtr<polyMesh> mesh = reader.mesh(runTime);
    reader.writeMesh(mesh(), format);


    if (args.optionFound("check"))
    {
        checkFireEdges(mesh());
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}

// ************************************************************************* //

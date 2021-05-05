/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

Description
    Code chunk for converting face and point sets - included by foamToVTK.

\*---------------------------------------------------------------------------*/

// flag to top-level code to signal early stop.
bool wroteTopoSet = false;

// Write faceSet (as PolyData)
if (faceSetName.size())
{
    // Load
    faceSet set(mesh, faceSetName);

    fileName outputName
    (
        outputDir/regionDir/"face-set"
      / set.name()/set.name() + timeDesc
    );

    Info<< "    faceSet   : "
        << args.relativePath(outputName) << nl;

    vtk::writeFaceSet
    (
        mesh,
        set,
        writeOpts,
        outputName,
        Pstream::parRun()
    );

    wroteTopoSet = true;
}


// Write pointSet (as PolyData)
if (pointSetName.size())
{
    // Load
    pointSet set(mesh, pointSetName);

    fileName outputName
    (
        outputDir/regionDir/"point-set"
      / set.name()/set.name() + timeDesc
    );

    Info<< "    pointSet  : "
        << args.relativePath(outputName) << nl;

    vtk::writePointSet
    (
        mesh,
        set,
        writeOpts,
        outputName,
        Pstream::parRun()
    );

    wroteTopoSet = true;
}


// ************************************************************************* //

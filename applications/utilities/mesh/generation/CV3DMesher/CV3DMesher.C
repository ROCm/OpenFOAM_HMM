/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2007 OpenCFD Ltd.
    \\/      M anipulation   |
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
    CV3DMesher

Description
    Conformal-Voronoi 3D automatic mesher

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CV3D.H"
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("surface");
    argList::validOptions.insert("pointsFile", "<filename>");

#   include "setRootCase.H"
#   include "createTime.H"

    // Read control dictionary
    // ~~~~~~~~~~~~~~~~~~~~~~~

    dictionary controlDict
    (
        IFstream(runTime.system()/args.executable() + "Dict")()
    );

    label nIterations(readLabel(controlDict.lookup("nIterations")));

    // Read the surface to conform to
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    querySurface surf
    (
        args.args()[1],
        readScalar(controlDict.lookup("includedAngle"))
    );

    Info<< nl
        << "Read surface with " << surf.size() << " triangles from file "
        << args.args()[1] << nl << endl;

    // Read and triangulation
    // ~~~~~~~~~~~~~~~~~~~~~~

    CV3D mesh(controlDict, surf);

    if (args.options().found("pointsFile"))
    {
        fileName pointFileName(args.options()["pointsFile"]);
        mesh.insertPoints(pointFileName);
        mesh.insertSurfacePointPairs();
        mesh.boundaryConform();
    }
    else
    {
        mesh.insertGrid();
        mesh.insertSurfacePointPairs();
        mesh.boundaryConform();
    }

    for (int iter=1; iter<=nIterations; iter++)
    {
        Info<< nl
            << "Relaxation iteration " << iter << nl
            << "~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

        mesh.removeSurfacePointPairs();
        mesh.insertSurfacePointPairs();
        mesh.boundaryConform();
    }

    mesh.write();

    mesh.writeDual("dualMesh.obj");

    mesh.writeMesh(runTime);

    Info << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

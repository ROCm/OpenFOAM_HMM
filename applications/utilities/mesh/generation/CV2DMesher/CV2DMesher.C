/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2007-2009 OpenCFD Ltd.
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
    CV2DMesher

Description
    Conformal-Voronoi 2D automatic mesher with grid or read initial points
    and point position relaxation with optional "squarification".

\*---------------------------------------------------------------------------*/

#include "CV2D.H"
#include "argList.H"
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

    argList args(argc, argv);

    // Read control dictionary
    // ~~~~~~~~~~~~~~~~~~~~~~~
    dictionary controlDict(IFstream(args.executable() + "Dict")());

    label nIterations(readLabel(controlDict.lookup("nIterations")));


    // Read the surface to conform to
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    querySurface surf(args.args()[1]);
    surf.writeTreeOBJ();

    Info<< nl
        << "Read surface with " << surf.size() << " triangles from file "
        << args.args()[1] << nl << endl;


    // Read and triangulation
    // ~~~~~~~~~~~~~~~~~~~~~~

    CV2D mesh(controlDict, surf);
    if (args.options().found("pointsFile"))
    {
        fileName pointFileName(IStringStream(args.options()["pointsFile"])());
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

        scalar relax =
           mesh.meshingControls().relaxationFactorStart
         +
         (
            mesh.meshingControls().relaxationFactorEnd
          - mesh.meshingControls().relaxationFactorStart
         )
         *scalar(iter)/scalar(nIterations);

        mesh.newPoints(relax);
        mesh.removeSurfacePointPairs();
        mesh.insertSurfacePointPairs();
        mesh.boundaryConform();
    }

    mesh.write();

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

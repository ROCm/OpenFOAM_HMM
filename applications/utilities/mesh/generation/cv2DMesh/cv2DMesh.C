/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2007-2010 OpenCFD Ltd.
    \\/      M anipulation   |
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
    cv2DMesh

Description
    Conformal-Voronoi 2D extruding automatic mesher with grid or read
    initial points and point position relaxation with optional
    "squarification".

\*---------------------------------------------------------------------------*/

#include "CV2D.H"
#include "argList.H"
#include "IFstream.H"

#include "MeshedSurfaces.H"
#include "shortEdgeFilter2D.H"
#include "extrude2DMesh.H"
#include "polyMesh.H"
#include "PatchTools.H"
#include "patchTo2DpolyMesh.H"
#include "extrudeModel.H"
#include "polyTopoChange.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("surface");
    argList::validOptions.insert("pointsFile", "<filename>");

    #include "setRootCase.H"
    #include "createTime.H"

    // Read control dictionary
    // ~~~~~~~~~~~~~~~~~~~~~~~
    dictionary controlDict(IFstream("system/" + args.executable() + "Dict")());
    dictionary shortEdgeFilterDict(controlDict.subDict("shortEdgeFilter"));
    dictionary extrusionDict(controlDict.subDict("extrusion"));

    Switch extrude(extrusionDict.lookup("extrude"));
    label nIterations(readLabel(controlDict.lookup("nIterations")));
    label sefDebug(shortEdgeFilterDict.lookupOrDefault<label>("debug", 0));

    // Read the surface to conform to
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    querySurface surf(args.args()[1]);
//    surf.writeTreeOBJ();

//    Info<< nl
//        << "Read surface with " << surf.size() << " triangles from file "
//        << args.args()[1] << nl << endl;

   // surf.write("surface.obj");

    // Read and triangulation
    // ~~~~~~~~~~~~~~~~~~~~~~
    CV2D mesh(runTime, controlDict);
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

        Info<< "Relaxation = " << relax << endl;

        mesh.newPoints(relax);
    }

    mesh.write();

    Info<< "Finished Delaunay in = "
        << runTime.cpuTimeIncrement() << " s." << endl;

    shortEdgeFilter2D sef(mesh, shortEdgeFilterDict);

    shortEdgeFilter2D::debug = sefDebug;

    sef.filter();

    Info<< "Meshed surface after edge filtering :" << endl;
    sef.fMesh().writeStats(Info);

    Info<< "Write .obj file : MeshedSurface.obj" << endl;
    sef.fMesh().write("MeshedSurface.obj");

    Info<< "Finished filtering in = "
        << runTime.cpuTimeIncrement() << " s." << endl;

    patchTo2DpolyMesh poly2DMesh
    (
        sef.fMesh(),
        sef.patchNames(),
        sef.patchSizes(),
        sef.mapEdgesRegion()
    );

    poly2DMesh.createMesh();

    polyMesh pMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        xferMove(poly2DMesh.points()),
        xferMove(poly2DMesh.faces()),
        xferMove(poly2DMesh.owner()),
        xferMove(poly2DMesh.neighbour())
    );

    Info<< "Constructing patches." << endl;
    List<polyPatch*> patches(poly2DMesh.patchNames().size());

    forAll(patches, patchI)
    {
        patches[patchI] = new polyPatch
        (
            poly2DMesh.patchNames()[patchI],
            poly2DMesh.patchSizes()[patchI],
            poly2DMesh.patchStarts()[patchI],
            patchI,
            pMesh.boundaryMesh()
        );
    }

    pMesh.addPatches(patches);

    if (extrude)
    {
        // Point generator
        autoPtr<extrudeModel> model(extrudeModel::New(extrusionDict));

        extrude2DMesh extruder(pMesh, extrusionDict, model());

        polyTopoChange meshMod(pMesh.boundaryMesh().size());

        extruder.addFrontBackPatches();
        extruder.setRefinement(meshMod);

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(pMesh, false);

        pMesh.updateMesh(morphMap);

        pMesh.setInstance("constant");
    }

    pMesh.write();

    Info<< "Finished extruding in = "
        << runTime.cpuTimeIncrement() << " s." << endl;

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

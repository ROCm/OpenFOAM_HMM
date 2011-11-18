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

#include "MeshedSurfaces.H"
#include "shortEdgeFilter2D.H"
#include "extrude2DMesh.H"
#include "polyMesh.H"
#include "patchToPoly2DMesh.H"
#include "extrudeModel.H"
#include "polyTopoChange.H"
#include "edgeCollapser.H"
#include "relaxationModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validOptions.insert("pointsFile", "<filename>");

    #include "addOverwriteOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    // Read control dictionary
    // ~~~~~~~~~~~~~~~~~~~~~~~
    dictionary controlDict(IFstream("system/" + args.executable() + "Dict")());
    dictionary shortEdgeFilterDict(controlDict.subDict("shortEdgeFilter"));
    dictionary extrusionDict(controlDict.subDict("extrusion"));

    Switch extrude(extrusionDict.lookup("extrude"));
    const bool overwrite = args.optionFound("overwrite");

    autoPtr<relaxationModel> relax
    (
        relaxationModel::New
        (
            controlDict.subDict("motionControl"),
            runTime
        )
    );

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

    while (runTime.loop())
    {
        Info<< nl << "Time = " << runTime.timeName() << endl;

        Info<< "Relaxation = " << relax->relaxation() << endl;

        mesh.newPoints(relax->relaxation());
    }

    mesh.write();

    Info<< "Finished Delaunay in = "
        << runTime.cpuTimeIncrement() << " s." << endl;

    Info<< "Begin filtering short edges:" << endl;
    shortEdgeFilter2D sef(mesh, shortEdgeFilterDict);

    sef.filter();

    Info<< "Meshed surface after edge filtering :" << endl;
    sef.fMesh().writeStats(Info);

    Info<< "Write .obj file of the 2D mesh: MeshedSurface.obj" << endl;
    sef.fMesh().write("MeshedSurface.obj");

    Info<< "Finished filtering in = "
        << runTime.cpuTimeIncrement() << " s." << endl;

    Info<< "Begin constructing a polyMesh:" << endl;

    patchToPoly2DMesh poly2DMesh
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
        Info<< "Begin extruding the polyMesh:" << endl;

        {
            // Point generator
            autoPtr<extrudeModel> model(extrudeModel::New(extrusionDict));

            extrude2DMesh extruder(pMesh, extrusionDict, model());

            extruder.addFrontBackPatches();

            polyTopoChange meshMod(pMesh.boundaryMesh().size());

            extruder.setRefinement(meshMod);

            autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(pMesh, false);

            pMesh.updateMesh(morphMap);
        }

        {
            edgeCollapser collapser(pMesh);

            const edgeList& edges = pMesh.edges();
            const pointField& points = pMesh.points();

            const boundBox& bb = pMesh.bounds();
            const scalar mergeDim = 1E-4 * bb.minDim();

            forAll(edges, edgeI)
            {
                const edge& e = edges[edgeI];

                scalar d = e.mag(points);

                if (d < mergeDim)
                {
                    Info<< "Merging edge " << e << " since length " << d
                        << " << " << mergeDim << endl;

                    // Collapse edge to e[0]
                    collapser.collapseEdge(edgeI, e[0]);
                }
            }

            polyTopoChange meshModCollapse(pMesh);

            collapser.setRefinement(meshModCollapse);

            // Create a mesh from topo changes.
            autoPtr<mapPolyMesh> morphMap =
                meshModCollapse.changeMesh(pMesh, false);

            pMesh.updateMesh(morphMap);
        }

        if (!overwrite)
        {
            runTime++;
        }
        else
        {
            pMesh.setInstance("constant");
        }

    }

    pMesh.write();

    Info<< "Finished extruding in = "
        << runTime.cpuTimeIncrement() << " s." << endl;

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

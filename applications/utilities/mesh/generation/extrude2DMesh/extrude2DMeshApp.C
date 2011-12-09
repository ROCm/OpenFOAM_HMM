/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    extrude2DMesh

Description
    Takes 2D mesh (all faces 2 points only, no front and back faces) and
    creates a 3D mesh by extruding with specified thickness.

Note
    Not sure about the walking of the faces to create the front and back faces.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "extrude2DMesh.H"
#include "extrudeModel.H"
#include "polyTopoChange.H"
#include "MeshedSurface.H"
#include "edgeCollapser.H"
#include "addPatchCellLayer.H"
#include "patchToPoly2DMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

enum ExtrudeMode
{
    POLYMESH2D,
    MESHEDSURFACE
};

namespace Foam
{
    template<>
    const char* NamedEnum<ExtrudeMode, 2>::names[] =
    {
        "polyMesh2D",
        "MeshedSurface"
    };
}

static const NamedEnum<ExtrudeMode, 2> ExtrudeModeNames;


//pointField moveInitialPoints
//(
//    primitiveFacePatch& fMesh,
//    const extrudeModel& model
//)
//{
//    pointField layer0Points(fMesh.nPoints());
//    pointField layer1Points(fMesh.nPoints());
//    pointField displacement(fMesh.nPoints());

//    forAll(layer0Points, pointI)
//    {
//        const labelList& meshPoints = fMesh.meshPoints();
//        label meshPointI = meshPoints[pointI];

//        layer0Points[meshPointI] = model
//        (
//            fMesh.points()[meshPointI],
//            fMesh.pointNormals()[pointI],
//            0
//        );

//        layer1Points[meshPointI] = model
//        (
//            fMesh.points()[meshPointI],
//            fMesh.pointNormals()[pointI],
//            1
//        );

//        displacement[pointI] =
//            layer1Points[meshPointI]
//          - layer0Points[meshPointI];
//    }

//    fMesh.movePoints(layer0Points);

//    return displacement;
//}


// Main program:

int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"

    argList::validArgs.append("surfaceFormat");

    #include "setRootCase.H"

    Info<< "Create time\n" << endl;

    Time runTimeExtruded
    (
        Time::controlDictName,
        args.rootPath(),
        args.caseName()
    );

    runTimeExtruded.functionObjects().off();

    const bool overwrite = args.optionFound("overwrite");

    const ExtrudeMode surfaceFormat = ExtrudeModeNames[args[1]];

    Info<< "Extruding from " << ExtrudeModeNames[surfaceFormat] << endl;

    IOdictionary extrude2DMeshDict
    (
        IOobject
        (
            "extrude2DMeshDict",
            runTimeExtruded.system(),
            runTimeExtruded,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // Point generator
    autoPtr<extrudeModel> model(extrudeModel::New(extrude2DMeshDict));

    autoPtr<MeshedSurface<face> > fMesh;

    autoPtr<polyMesh> mesh;

    autoPtr<polyTopoChange> meshMod;

    labelListList extrudeEdgePatches;

    if (surfaceFormat == MESHEDSURFACE)
    {
        fMesh.set(new MeshedSurface<face>("MeshedSurface.obj"));

        EdgeMap<label> edgeRegionMap;
        wordList patchNames(1, "default");
        labelList patchSizes(1, fMesh->nEdges() - fMesh->nInternalEdges());

        const edgeList& edges = fMesh->edges();
        forAll(edges, edgeI)
        {
            if (!fMesh->isInternalEdge(edgeI))
            {
                edgeRegionMap.insert(edges[edgeI], 0);
            }
        }

        patchToPoly2DMesh poly2DMesh
        (
            fMesh(),
            patchNames,
            patchSizes,
            edgeRegionMap
        );

        poly2DMesh.createMesh();

        mesh.set
        (
            new polyMesh
            (
                IOobject
                (
                    polyMesh::defaultRegion,
                    runTimeExtruded.constant(),
                    runTimeExtruded,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                xferMove(poly2DMesh.points()),
                xferMove(poly2DMesh.faces()),
                xferMove(poly2DMesh.owner()),
                xferMove(poly2DMesh.neighbour())
            )
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
                mesh->boundaryMesh()
            );
        }

        mesh->addPatches(patches);
    }
    else if (surfaceFormat == POLYMESH2D)
    {
        mesh.set
        (
            new polyMesh
            (
                Foam::IOobject
                (
                    Foam::polyMesh::defaultRegion,
                    runTimeExtruded.timeName(),
                    runTimeExtruded,
                    Foam::IOobject::MUST_READ
                )
            )
        );
    }

    // Engine to extrude mesh
    extrude2DMesh extruder(mesh(), extrude2DMeshDict, model());

    extruder.addFrontBackPatches();

    meshMod.set(new polyTopoChange(mesh().boundaryMesh().size()));

    extruder.setRefinement(meshMod());

    // Create a mesh from topo changes.
    autoPtr<mapPolyMesh> morphMap = meshMod->changeMesh(mesh(), false);

    mesh->updateMesh(morphMap);

    {
        edgeCollapser collapser(mesh());

        const edgeList& edges = mesh->edges();
        const pointField& points = mesh->points();

        const boundBox& bb = mesh->bounds();
        const scalar mergeDim = 1E-4 * bb.minDim();

        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            scalar d = e.mag(points);

            if (d < mergeDim)
            {
                Info<< "Merging edge " << e << " since length " << d
                    << " << " << mergeDim << nl;

                // Collapse edge to e[0]
                collapser.collapseEdge(edgeI, e[0]);
            }
        }

        polyTopoChange meshModCollapse(mesh());

        collapser.setRefinement(meshModCollapse);

        // Create a mesh from topo changes.
        autoPtr<mapPolyMesh> morphMap
            = meshModCollapse.changeMesh(mesh(), false);

        mesh->updateMesh(morphMap);
    }


    mesh->setInstance(runTimeExtruded.constant());

    // Take over refinement levels and write to new time directory.
    Pout<< "Writing extruded mesh to time = " << runTimeExtruded.timeName()
        << nl << endl;

    mesh().write();

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

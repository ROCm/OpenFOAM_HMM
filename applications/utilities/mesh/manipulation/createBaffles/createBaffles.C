/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Description
    Makes internal faces into boundary faces. Does not duplicate points. Use
    mergeOrSplitBaffles if you want this.

    Note: if any coupled patch face is selected for baffling automatically
    the opposite member is selected for baffling as well. Note that this
    is the same as repatching. This was added only for convenience so
    you don't have to filter coupled boundary out of your set.

\*---------------------------------------------------------------------------*/

#include "syncTools.H"
#include "argList.H"
#include "Time.H"
#include "faceSet.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"
#include "polyAddFace.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("set");
    argList::validArgs.append("patch");
    argList::validOptions.insert("additionalPatches", "(patch2 .. patchN)");
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const faceZoneMesh& faceZones = mesh.faceZones();

    // Faces to baffle
    word setName(args.additionalArgs()[0]);
    Info<< "Reading faceSet from " << setName << nl << endl;
    faceSet facesToSplit(mesh, setName);
    // Make sure set is synchronised across couples
    facesToSplit.sync(mesh);
    Info<< "Read " << returnReduce(facesToSplit.size(), sumOp<label>())
        << " faces from " << setName << nl << endl;


    // Patches to put baffles into
    labelList newPatches(1);

    word patchName(args.additionalArgs()[1]);
    newPatches[0] = patches.findPatchID(patchName);
    Info<< "Using patch " << patchName
        << " at index " << newPatches[0] << endl;

    if (newPatches[0] == -1)
    {
        FatalErrorIn(args.executable())
            << "Cannot find patch " << patchName << endl
            << "Valid patches are " << patches.names() << exit(FatalError);
    }


    // Additional patches
    if (args.options().found("additionalPatches"))
    {
        const wordList patchNames
        (
            IStringStream(args.options()["additionalPatches"])()
        );

        forAll(patchNames, i)
        {
            label patchI = patches.findPatchID(patchNames[i]);
            Info<< "Using additional patch " << patchNames[i]
                << " at index " << patchI << endl;

            if (patchI == -1)
            {
                FatalErrorIn(args.executable())
                    << "Cannot find patch " << patchNames[i] << endl
                    << "Valid patches are " << patches.names()
                    << exit(FatalError);
            }

            newPatches.append(patchI);
        }
    }


    bool overwrite = args.options().found("overwrite");



    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);


    // Mesh change container
    polyTopoChange meshMod(mesh);


    // Do the actual changes
    // Note order in which faces are modified/added is so they end up correctly
    // for cyclic patches (original faces first and then reversed faces)
    // since otherwise it will have trouble matching baffles.

    label nBaffled = 0;
 
    forAll(newPatches, i)
    {
        label newPatchI = newPatches[i];

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            if (facesToSplit.found(faceI))
            {
                const face& f = mesh.faces()[faceI];
                label zoneID = faceZones.whichZone(faceI);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                }

                if (i == 0)
                {
                    // First usage of face. Modify.
                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            f,                          // modified face
                            faceI,                      // label of face
                            mesh.faceOwner()[faceI],    // owner
                            -1,                         // neighbour
                            false,                      // face flip
                            newPatchI,                  // patch for face
                            false,                      // remove from zone
                            zoneID,                     // zone for face
                            zoneFlip                    // face flip in zone
                        )
                    );
                }
                else
                {
                    // Second or more usage of face. Add.
                    meshMod.setAction
                    (
                        polyAddFace
                        (
                            f,                          // modified face
                            mesh.faceOwner()[faceI],    // owner
                            -1,                         // neighbour
                            -1,                         // master point
                            -1,                         // master edge
                            faceI,                      // master face
                            false,                      // face flip
                            newPatchI,                  // patch for face
                            zoneID,                     // zone for face
                            zoneFlip                    // face flip in zone
                        )
                    );
                }

                nBaffled++;
            }
        }

        // Add the reversed face.
        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            if (facesToSplit.found(faceI))
            {
                const face& f = mesh.faces()[faceI];
                label zoneID = faceZones.whichZone(faceI);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                }
                label nei = mesh.faceNeighbour()[faceI];

                meshMod.setAction
                (
                    polyAddFace
                    (
                        f.reverseFace(),            // modified face
                        nei,                        // owner
                        -1,                         // neighbour
                        -1,                         // masterPointID
                        -1,                         // masterEdgeID
                        faceI,                      // masterFaceID,
                        true,                       // face flip
                        newPatchI,                  // patch for face
                        zoneID,                     // zone for face
                        zoneFlip                    // face flip in zone
                    )
                );

                nBaffled++;
            }
        }

        // Modify any boundary faces
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (patches[newPatchI].coupled() && pp.coupled())
            {
                // Do not allow coupled faces to be moved to different coupled
                // patches.
            }
            else
            {
                forAll(pp, i)
                {
                    label faceI = pp.start()+i;

                    if (facesToSplit.found(faceI))
                    {
                        const face& f = mesh.faces()[faceI];
                        label zoneID = faceZones.whichZone(faceI);
                        bool zoneFlip = false;
                        if (zoneID >= 0)
                        {
                            const faceZone& fZone = faceZones[zoneID];
                            zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                        }

                        if (i == 0)
                        {
                            // First usage of face. Modify.
                            meshMod.setAction
                            (
                                polyModifyFace
                                (
                                    f,                      // modified face
                                    faceI,                  // label of face
                                    mesh.faceOwner()[faceI],// owner
                                    -1,                     // neighbour
                                    false,                  // face flip
                                    newPatchI,              // patch for face
                                    false,                  // remove from zone
                                    zoneID,                 // zone for face
                                    zoneFlip                // face flip in zone
                                )
                            );
                        }
                        else
                        {
                            // Second or more usage of face. Add.
                            meshMod.setAction
                            (
                                polyAddFace
                                (
                                    f,                      // modified face
                                    mesh.faceOwner()[faceI],// owner
                                    -1,                     // neighbour
                                    -1,                     // master point
                                    -1,                     // master edge
                                    faceI,                  // master face
                                    false,                  // face flip
                                    newPatchI,              // patch for face
                                    zoneID,                 // zone for face
                                    zoneFlip                // face flip in zone
                                )
                            );
                        }

                        nBaffled++;
                    }
                }
            }
        }
    }


    Info<< "Converted " << returnReduce(nBaffled, sumOp<label>())
        << " faces into boundary faces on patch " << patchName << nl << endl;

    if (!overwrite)
    {
        runTime++;
    }

    // Change the mesh. Change points directly (no inflation).
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    Info<< "Writing mesh to " << runTime.timeName() << endl;

    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

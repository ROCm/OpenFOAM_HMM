/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Utility to create patches out of selected boundary faces. Faces come either
    from existing patches or from a faceSet.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "SortableList.H"
#include "OFstream.H"
#include "meshTools.H"
#include "faceSet.H"
#include "IOPtrList.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<dictionary>, 0);
}


label getPatch(const polyBoundaryMesh& patches, const word& patchName)
{
    label patchI = patches.findPatchID(patchName);

    if (patchI == -1)
    {
        FatalErrorIn("createPatch")
            << "Cannot find source patch " << patchName
            << endl << "Valid patch names are " << patches.names()
            << exit(FatalError);
    }

    return patchI;
}


void changePatchID
(
    const polyMesh& mesh,
    const label faceID,
    const label patchID,
    polyTopoChange& meshMod
)
{
    const label zoneID = mesh.faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            mesh.faces()[faceID],               // face
            faceID,                             // face ID
            mesh.faceOwner()[faceID],           // owner
            -1,                                 // neighbour
            false,                              // flip flux
            patchID,                            // patch ID
            false,                              // remove from zone
            zoneID,                             // zone ID
            zoneFlip                            // zone flip
        )
    );
}


// Filter out the empty patches.
void filterPatches(polyMesh& mesh)
{
    Info<< "Filtering empty patches." << endl;

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Old patches.
    DynamicList<polyPatch*> allPatches(patches.size());

    labelList compactPatchMap(patches.size());

    // Copy old patches.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.size() > 0)
        {
            compactPatchMap[patchI] = allPatches.size();

            allPatches.append
            (
                pp.clone
                (
                    patches,
                    allPatches.size(),
                    pp.size(),
                    pp.start()
                ).ptr()
            );
        }
        else
        {
            Info<< "Removing empty patch " << pp.name() << " at position "
                << patchI << endl;

            compactPatchMap[patchI] = -1;
        }
    }
    allPatches.shrink();

    mesh.removeBoundary();
    mesh.addPatches(allPatches);
}


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"


    Info<< "Reading createPatchDict\n" << endl;

    PtrList<dictionary> patchSources
    (
        IOdictionary
        (
            IOobject
            (
                "createPatchDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("patches")
    );

    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // 1. All all new patches
    // ~~~~~~~~~~~~~~~~~~~~~~

    // Old and new patches.
    DynamicList<polyPatch*> allPatches(patches.size() + patchSources.size());

    // Copy old patches.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        Info<< "Copying patch " << pp.name() << " at position "
            << patchI << endl;

        allPatches.append
        (
            pp.clone
            (
                patches,
                patchI,
                pp.size(),
                pp.start()
            ).ptr()
        );
    }

    forAll(patchSources, addedI)
    {
        const dictionary& dict = patchSources[addedI];

        word patchName(dict.lookup("name"));

        label destPatchI = patches.findPatchID(patchName);

        word patchType(dict.lookup("type"));

        if (destPatchI == -1)
        {
            destPatchI = allPatches.size();

            Info<< "Adding new patch " << patchName << " of type " << patchType
                << " as patch " << destPatchI << endl;

            // Add an empty patch.
            allPatches.append
            (
                polyPatch::New
                (
                    patchType,
                    patchName,
                    0,
                    mesh.nFaces(),
                    destPatchI,
                    patches
                ).ptr()
            );
        }
    }

    allPatches.shrink();

    mesh.removeBoundary();
    mesh.addPatches(allPatches);



    // 2. Repatch faces
    // ~~~~~~~~~~~~~~~~

    polyTopoChange meshMod(mesh);


    forAll(patchSources, addedI)
    {
        const dictionary& dict = patchSources[addedI];

        word patchName(dict.lookup("name"));

        label destPatchI = patches.findPatchID(patchName);

        if (destPatchI == -1)
        {
            FatalErrorIn(args.executable()) << "patch " << patchName
                << " not added. Problem." << abort(FatalError);
        }

        word sourceType(dict.lookup("constructFrom"));

        if (sourceType == "patches")
        {
            wordList patchSources(dict.lookup("patches"));

            // Repatch faces of the patches.
            forAll(patchSources, sourceI)
            {
                label patchI = getPatch(patches, patchSources[sourceI]);

                const polyPatch& pp = patches[patchI];

                Info<< "Moving faces from patch " << pp.name()
                    << " to patch " << destPatchI << endl;

                forAll(pp, i)
                {
                    changePatchID
                    (
                        mesh,
                        pp.start() + i,
                        destPatchI,
                        meshMod
                    );
                }
            }
        }
        else if (sourceType == "set")
        {
            word setName(dict.lookup("set"));

            faceSet faces(mesh, setName);

            Info<< "Read " << faces.size() << " faces from faceSet "
                << faces.name() << endl;

            // Sort (since faceSet contains faces in arbitrary order)
            labelList faceLabels(faces.toc());

            SortableList<label> patchFaces(faceLabels);

            forAll(patchFaces, i)
            {
                label faceI = patchFaces[i];

                if (mesh.isInternalFace(faceI))
                {
                    FatalErrorIn(args.executable())
                        << "Face " << faceI << " specified in set "
                        << faces.name()
                        << " is not an external face of the mesh." << endl
                        << "This application can only repatch existing boundary"
                        << " faces." << exit(FatalError);
                }

                changePatchID
                (
                    mesh,
                    faceI,
                    destPatchI,
                    meshMod
                );
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Invalid source type " << sourceType << endl
                << "Valid source types are 'patches' 'set'" << exit(FatalError);
        }
    }

    // Change mesh, no inflation
    meshMod.changeMesh(mesh, false);



    // 3. Remove zeros-sized patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    filterPatches(mesh);



    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    runTime++;

    // Write resulting mesh
    Info<< "Writing repatched mesh to " << runTime.timeName() << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

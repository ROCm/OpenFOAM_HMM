/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Utility to reorder cyclic and processor patches.

    Uses dummy morph to sort things out.

    Is bit of hack since polyMesh constructor already checks for coupled
    face areas so might bomb out. You might need to compile in a local
    processorPolyPatch that gives a warning instead to make this utility
    complete.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "OFstream.H"
#include "coupledPolyPatch.H"
#include "PstreamReduceOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    Info<< "Using geometry to calculate face correspondence across"
        << " coupled boundaries (processor, cyclic)" << nl
        << "This will only work for cyclics if they are parallel or"
        << " their rotation is defined across the origin" << nl
        << endl;

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    bool hasCoupled = false;
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            hasCoupled = true;

            break;
        }
    }

    reduce(hasCoupled, orOp<bool>());

    if (hasCoupled)
    {
        Info<< "Mesh has coupled patches ..." << nl << endl;

        // Dummy topo changes container
        polyTopoChange meshMod(mesh);

        // Do all changes
        Info<< "Doing dummy mesh morph to correct face ordering ..."
            << endl;

        runTime++;

        faceList oldFaces(mesh.faces());

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(morphMap);

        // Move mesh (since morphing does not do this)
        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }

        // Find out if anything changed.
        // Problem is that - if the topo change is only rotation -
        // face::operator= does not detect a change so use labelList::operator=
        bool meshChanged = false;

        if (mesh.faces().size() != oldFaces.size())
        {
            meshChanged = true;
        }
        else
        {
            forAll(mesh.faces(), faceI)
            {
                const labelList& f = mesh.faces()[faceI];
                const labelList& oldF = oldFaces[faceI];

                if (f != oldF)
                {
                    meshChanged = true;
                    break;
                }
            }
        }

        reduce(meshChanged, orOp<label>());

        if (meshChanged)
        {
            // Set the precision of the points data to 10
            IOstream::defaultPrecision(10);

            // Write resulting mesh
            Info << "Writing morphed mesh to time " << runTime.value() << endl;

            mesh.write();
        }
        else
        {
            Info << "Mesh ordering ok. Nothing changed." << endl;
        }
    }
    else
    {
        Info<< "Mesh has no coupled patches. Nothing changed ..." << nl << endl;
    }
        

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //

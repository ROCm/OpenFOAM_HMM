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

Application
    snappyHexMesh

Description
    Automatic split hex mesher. Refines and snaps to surface.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "autoHexMeshDriver.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createMesh.H"

    Info<< "Read mesh in = "
        << runTime.cpuTimeIncrement() << " s" << endl;

    // Read decomposePar dictionary
    IOdictionary decomposeDict
    (
        IOobject
        (
            "decomposeParDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read meshing dictionary
    IOdictionary meshDict
    (
       IOobject
       (
            "snappyHexMeshDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
       )
    );

    // refinement parameters
    const dictionary& refineDict = meshDict.subDict("refineDict");

    // snap-to-surface parameters
    const dictionary& snapDict = meshDict.subDict("snapDict");

    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("motionDict");

    // layer addition parameters
    const dictionary& layerDict = meshDict.subDict("layerDict");


    // Main meshing driver. Read surfaces. Determine initial intersections.
    autoHexMeshDriver meshEngine
    (
        mesh,
        meshDict,       // global control parameters
        refineDict,     // refinement parameters
        decomposeDict
    );

    Switch wantRefine(meshDict.lookup("doRefine"));
    Switch wantSnap(meshDict.lookup("doSnap"));
    Switch wantLayers(meshDict.lookup("doLayers"));

    if (wantRefine)
    {
        meshEngine.doRefine(refineDict, wantSnap);
    }

    if (wantSnap)
    {
        meshEngine.doSnap(snapDict, motionDict);
    }

    if (wantLayers)
    {
        meshEngine.doLayers(layerDict, motionDict);
    }

    Info<< "Finished meshing in = "
        << runTime.elapsedCpuTime() << " s." << endl;

    Pout<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

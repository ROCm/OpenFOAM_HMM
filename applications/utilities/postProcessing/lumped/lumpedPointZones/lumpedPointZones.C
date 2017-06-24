/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    lumpedPointZones

Description
    Produces a VTK PolyData file \c lumpedPointZones.vtp in which the
    segmentation of the pressure integration zones can be visualized
    for diagnostic purposes. Does not use external coupling.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"

#include "lumpedPointTools.H"
#include "lumpedPointIOMovement.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create lumpedPointZones.vtp to verify the segmentation of "
        "pressure integration zones used by lumpedPoint BC."
    );
    argList::noParallel();    // The VTP writer is not yet in parallel

    argList::noFunctionObjects();  // Never use function objects
    argList::addBoolOption
    (
        "verbose",
        "increased verbosity"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"

    // const bool verbose = args.optionFound("verbose");

    #include "createTime.H"

    runTime.setTime(instant(0, runTime.constant()), 0);

    #include "createNamedPolyMesh.H"

    autoPtr<lumpedPointIOMovement> movement = lumpedPointIOMovement::New
    (
        runTime
    );

    if (!movement.valid())
    {
        Info<< "no valid movement found" << endl;
        return 1;
    }

    const labelList patchLst = lumpedPointTools::lumpedPointPatchList(mesh);
    if (patchLst.empty())
    {
        Info<< "no patch list found" << endl;
        return 2;
    }

    pointIOField points0 = lumpedPointTools::points0Field(mesh);
    movement().setMapping(mesh, patchLst, points0);

    // Initial geometry, but with zone colouring
    movement().writeZonesVTP("lumpedPointZones.vtp", mesh, points0);

    // Initial positions/rotations
    movement().writeStateVTP("initialState.vtp");

    Info<< nl
        << "wrote 'lumpedPointZones.vtp'" << nl
        << "wrote 'initialState.vtp'" << nl
        << "End\n" << endl;

    return 0;
}

// ************************************************************************* //

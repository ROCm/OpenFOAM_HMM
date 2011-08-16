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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "IStringStream.H"
#include "octree.H"
#include "octreeDataCell.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("point (x y z)");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    label nReps = 100000;

    const point sample = args.argRead<point>(1);

    treeBoundBox meshBb(mesh.bounds());

    // Calculate typical cell related size to shift bb by.
    scalar typDim = meshBb.avgDim()/(2.0*Foam::cbrt(scalar(mesh.nCells())));

    treeBoundBox shiftedBb
    (
        meshBb.min(),
        meshBb.max() + vector(typDim, typDim, typDim)
    );

    Info<< "Mesh" << endl;
    Info<< "   bounding box           : " << meshBb << endl;
    Info<< "   bounding box (shifted) : " << shiftedBb << endl;
    Info<< "   typical dimension      : " << shiftedBb.typDim() << endl;

    Info<< "Initialised mesh in "
        << runTime.cpuTimeIncrement() << " s" << endl;

    // Wrap indices and mesh information into helper object
    octreeDataCell shapes(mesh);

    octree<octreeDataCell> oc
    (
        shiftedBb,  // overall bounding box
        shapes,     // all information needed to do checks on cells
        1,          // min. levels
        10.0,       // max. size of leaves
        10.0        // maximum ratio of cubes v.s. cells
    );

    for (label i = 0; i < nReps - 1 ; i++)
    {
        oc.find(sample);
    }

    Info<< "Point:" << sample << " is in shape "
        << oc.find(sample) << endl;

    oc.printStats(Info);

    Info<< "Found in octree " << nReps << " times in "
        << runTime.cpuTimeIncrement() << " s" << endl;

    indexedOctree<treeDataCell> ioc
    (
        treeDataCell(true, mesh),
        shiftedBb,
        8,      // maxLevel
        10,     // leafsize
        3.0     // duplicity
    );

    for (label i = 0; i < nReps - 1 ; i++)
    {
        ioc.findInside(sample);
    }

    Info<< "Point:" << sample << " is in shape "
        << ioc.findInside(sample)
        << ", where the possible cells were:" << nl
        << ioc.findIndices(sample)
        << endl;

    Info<< "Found in indexedOctree " << nReps << " times in "
        << runTime.cpuTimeIncrement() << " s" << endl;

    for (label i = 0; i < nReps - 1 ; i++)
    {
        mesh.findCell(sample);
    }

    Info<< "Point:" << sample << " is in cell  "
        << mesh.findCell(sample) << endl;

    Info<< "Found in mesh.findCell " << nReps << " times in "
        << runTime.cpuTimeIncrement() << " s" << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

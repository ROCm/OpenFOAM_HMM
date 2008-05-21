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
    Sample field data with a choice of interpolation schemes, sampling options
    and write formats.

    Keywords:

    setFormat: set output format, choice of
        xmgr
        jplot
        gnuplot
        raw

    surfaceFormat: surface output format, choice of
        null        : suppress output
        foamFile    : separate points, faces and values file
        dx          : DX scalar or vector format
        vtk         : VTK ascii format
        raw         : x y z value format for use with e.g. gnuplot 'splot'.
        stl         : ascii stl. Does not contain values!

    interpolationScheme: interpolation scheme, choice of
        cell          : use cell-centre value; constant over cells (default)
        cellPoint     : use cell-centre and vertex values
   	    cellPointFace : use cell-centre, vertex and face values.
          1] vertex values determined from neighbouring cell-centre values
          2] face values determined using the current face interpolation scheme
             for the field (linear, limitedLinear, etc.)

    fields: list of fields to sample

    sets: list of sets to sample, choice of
        uniform             evenly distributed points on line
        face                one point per face intersection
        midPoint            one point per cell, inbetween two face intersections
        midPointAndFace     combination of face and midPoint

        curve               specified points, not nessecary on line, uses
                            tracking
        cloud               specified points, uses findCell

        Option axis: how to write point coordinate. Choice of
          - x/y/z: x/y/z coordinate only
          - xyz: three columns
            (probably does not make sense for anything but raw)
          - distance: distance from start of sampling line (if uses line)
            or distance from first specified sampling point

        Type specific options:
            uniform, face, midPoint, midPointAndFace : start and end coordinate
            uniform: extra number of sampling points
            curve, cloud: list of coordinates

    surfaces: list of surfaces to sample, choice of
        plane : values on plane defined by point, normal.
        patch : values on patch.

    Runs in parallel.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOsampledSets.H"
#include "IOsampledSurfaces.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    IOsampledSets sSets(mesh, "sampleDict", true);
    IOsampledSurfaces sSurfaces(mesh, "sampleDict", true);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        // Handle geometry/topology changes
        polyMesh::readUpdateState state = mesh.readUpdate();

        sSets.readUpdate(state);
        sSurfaces.readUpdate(state);

        sSets.write();
        sSurfaces.write();

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

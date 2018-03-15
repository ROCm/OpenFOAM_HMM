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

Description
    Test bounding box behaviour

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "boundBox.H"
#include "treeBoundBox.H"
#include "cellModel.H"

using namespace Foam;

//- simple helper to create a cube
boundBox cube(scalar start, scalar width)
{
    return boundBox
    (
        point::uniform(start),
        point::uniform(start + width)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    Info<<"boundBox faces: " << boundBox::faces << nl
        <<"hex faces: " << cellModel::ref(cellModel::HEX).modelFaces() << nl
        <<"tree-bb faces: " << treeBoundBox::faces << nl
        <<"tree-bb edges: " << treeBoundBox::edges << endl;

    boundBox bb = boundBox::greatBox;
    Info<<"great box: " << bb << endl;

    // bb.clear();
    // Info<<"zero box: " << bb << endl;

    bb = boundBox::invertedBox;
    Info<<"invalid box: " << bb << endl;
    Info<< nl << endl;

    if (Pstream::parRun())
    {
        bb = cube(Pstream::myProcNo(), 1.1);
        Pout<<"box: " << bb << endl;

        bb.reduce();
        Pout<<"reduced: " << bb << endl;
    }
    else
    {
        bb = cube(0, 1);
        Info<<"starting box: " << bb << endl;

        point pt(Zero);
        bb.add(pt);
        Info<<"enclose point " << pt << " -> " << bb << endl;

        pt = point(0,1.5,0.5);
        bb.add(pt);
        Info<<"enclose point " << pt << " -> " << bb << endl;

        pt = point(5,2,-2);
        bb.add(pt);
        Info<<"enclose point " << pt << " -> " << bb << endl;

        // restart with same points
        bb = boundBox::invertedBox;
        bb.add(point(1,1,1));
        bb.add(point::zero);
        bb.add(point(0,1.5,0.5));
        bb.add(point(5,2,-2));

        Info<<"repeated " << bb << endl;

        boundBox box1 = cube(0, 1);
        boundBox box2 = cube(0, 0.75);
        boundBox box3 = cube(0.5, 1);
        boundBox box4 = cube(-1, 0.5);

        Info<<"union of " << box1 << " and " << box2 << " => ";

        box1.add(box2);
        Info<< box1 << endl;

        box1.add(box3);
        Info<<"union with " << box3 << " => " << box1 << endl;

        box1.add(box4);
        Info<<"union with " << box4 << " => " << box1 << endl;
    }

    return 0;
}


// ************************************************************************* //

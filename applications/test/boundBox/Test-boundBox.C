/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
#include "cellModeller.H"

using namespace Foam;

//- simple helper to create a cube
boundBox cube(scalar start, scalar width)
{
    return boundBox
    (
        point(start, start, start),
        point(start + width, start + width, start + width)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    // #include "createTime.H"
    // #include "createMesh.H"

    const cellModel& hex = *(cellModeller::lookup("hex"));

    Info<<"boundBox faces: " << boundBox::faces << endl;
    Info<<"hex faces: " << hex.modelFaces() << endl;

    boundBox bb = boundBox::greatBox;
    Info<<"great box: " << bb << endl;

    if (Pstream::parRun())
    {
        bb = cube(Pstream::myProcNo(), 1.1);
        Pout<<"box: " << bb << endl;

        bb.reduce();
        Pout<<"reduced: " << bb << endl;
    }

    return 0;
}


// ************************************************************************* //

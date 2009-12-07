/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2009-2009 OpenCFD Ltd.
    \\/      M anipulation   |
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
    momentOfInertiaTest

Description
    Calculates the inertia tensor and principal axes and moments of a
    test face.

\*---------------------------------------------------------------------------*/

#include "ListOps.H"
#include "face.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    label nPts = 6;

    pointField pts(nPts);

    pts[0] = point(4.495, 3.717, -4.112);
    pts[1] = point(4.421, 3.932, -4.112);
    pts[2] = point(4.379, 4.053, -4.112);
    pts[3] = point(4.301, 4.026, -4.300);
    pts[4] = point(4.294, 4.024, -4.317);
    pts[5] = point(4.409, 3.687, -4.317);

    scalar density = 1.0;

    face f(identity(nPts));

    point Cf = f.centre(pts);

    tensor J = tensor::zero;

    J = f.inertia(pts, Cf, density);

    vector eVal = eigenValues(J);

    tensor eVec = eigenVectors(J);

    Info<< nl << "Inertia tensor of test face " << J << nl
        << "eigenValues (principal moments) " << eVal << nl
        << "eigenVectors (principal axes) " << eVec
        << endl;

    OFstream str("momentOfInertiaTest.obj");

    Info<< nl << "Writing test face and scaled principal axes to "
        << str.name() << endl;

    forAll(pts, ptI)
    {
        meshTools::writeOBJ(str, pts[ptI]);
    }

    str << "l";

    forAll(f, fI)
    {
        str << ' ' << fI + 1;
    }

    str << " 1" << endl;

    scalar scale = mag(Cf - pts[f[0]])/eVal.component(findMin(eVal));

    meshTools::writeOBJ(str, Cf);
    meshTools::writeOBJ(str, Cf + scale*eVal.x()*eVec.x());
    meshTools::writeOBJ(str, Cf + scale*eVal.y()*eVec.y());
    meshTools::writeOBJ(str, Cf + scale*eVal.z()*eVec.z());

    for (label i = nPts + 1; i < nPts + 4; i++)
    {
        str << "l " << nPts + 1 << ' ' << i + 1 << endl;
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //

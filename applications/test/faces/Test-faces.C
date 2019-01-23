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

Application
    Test-faces

Description
    Simple tests for various faces

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "labelledTri.H"
#include "pointList.H"
#include "ListOps.H"

using namespace Foam;


template<class Face>
void faceInfo(const Face& f, const UList<point>& points)
{
    Info<< f
        << " points:" << f.points(points)
        << " normal:" << f.unitNormal(points);
}


template<class Face>
void testSign
(
    const Face& f,
    const UList<point>& points,
    const UList<point>& testPoints
)
{
    for (const point& p : testPoints)
    {
        Info<< "  point:" << p << " sign=" << f.sign(p, points) << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    pointList points1
    ({
        { 0, 0, 0},
        { -1, -1, 1},
        { 1, -1, -1},
        { 1, 1, -1},
        { -1, 1, 1}
    });

    pointList points2 = ListOps::create<point>
    (
        points1,
        [](const point& p){ return point(p.x(), p.y(), -p.z()); }
    );

    pointList testPoints
    ({
        { -2, -2, -2},
        { -2, -2, 2},
        { 0, 0, 0},
        { 2, 2, -2},
        { 2, 2, 2}
    });

    face f1{ 1, 2, 3, 4 };
    Info<< "face:"; faceInfo(f1, points1); Info << nl;
    testSign(f1, points1, testPoints);

    Info<< "face:"; faceInfo(f1, points2); Info << nl;
    testSign(f1, points2, testPoints);
    Info<< nl;

    triFace t1{ 1, 2, 3 };
    Info<< "triFace:"; faceInfo(t1, points1); Info << nl;
    testSign(t1, points1, testPoints);

    Info<< "triFace:"; faceInfo(t1, points2); Info << nl;
    testSign(t1, points2, testPoints);
    Info<< nl;


    f1 = t1;
    Info<< "face:" << f1 << nl;

    f1 = t1.triFaceFace();
    Info<< "face:" << f1 << nl;

    // expect these to fail
    const bool throwingError = FatalError.throwExceptions();
    try
    {
        labelledTri l1{ 1, 2, 3, 10, 24 };
        Info<< "labelled:" << l1 << nl;
    }
    catch (const Foam::error& err)
    {
        WarningInFunction
            << "Caught FatalError " << err << nl << endl;
    }
    FatalError.throwExceptions(throwingError);

    labelledTri l2{ 1, 2, 3 };
    Info<< "labelled:" << l2 << nl;

    labelledTri l3{ 1, 2, 3, 10 };
    Info<< "labelled:" << l3 << nl;

    t1.flip();
    l3.flip();

    Info<< "flip:" << t1 << nl;
    Info<< "flip:" << l3 << nl;

    return 0;
}


// ************************************************************************* //

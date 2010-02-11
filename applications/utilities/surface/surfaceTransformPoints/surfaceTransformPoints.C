/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Transform (scale/rotate) a surface.
    Like transformPoints but for surfaces.

    The rollPitchYaw option takes three angles (degrees):
    - roll (rotation about x) followed by
    - pitch (rotation about y) followed by
    - yaw (rotation about z)

    The yawPitchRoll does yaw followed by pitch followed by roll.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include "IFstream.H"
#include "boundBox.H"
#include "transformField.H"
#include "Pair.H"
#include "quaternion.H"
#include "mathematicalConstants.H"

#include "MeshedSurfaces.H"

using namespace Foam;
using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transform (scale/rotate) a surface. "
        "Like transformPoints but for surfaces."
    );
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");
    argList::addOption
    (
        "translate",
        "vector",
        "translate by the specified <vector> - eg, '(1 0 0)'"
    );
    argList::addOption
    (
        "rotate",
        "(vectorA vectorB)",
        "transform in terms of a rotation between <vectorA> and <vectorB> "
        "- eg, '( (1 0 0) (0 0 1) )'"
    );
    argList::addOption
    (
        "scale",
        "vector",
        "scale by the specified amount - eg, '(0.001 0.001 0.001)' for a "
        "uniform [mm] to [m] scaling"
    );
    argList::addOption
    (
        "rollPitchYaw",
        "vector",
        "transform in terms of '( roll pitch yaw )' in degrees"
    );
    argList::addOption
    (
        "yawPitchRoll",
        "vector",
        "transform in terms of '( yaw pitch roll )' in degrees"
    );
    argList args(argc, argv);

    fileName surfFileName(args.additionalArgs()[0]);

    Info<< "Reading surf from " << surfFileName << " ..." << endl;

    fileName outFileName(args.additionalArgs()[1]);

    Info<< "Writing surf to " << outFileName << " ..." << endl;


    if (args.options().empty())
    {
        FatalErrorIn(args.executable())
            << "No options supplied, please use one or more of "
               "-translate, -rotate or -scale options."
            << exit(FatalError);
    }

    meshedSurface surf1(surfFileName);

    pointField points(surf1.points());

    if (args.optionFound("translate"))
    {
        vector transVector(args.optionLookup("translate")());

        Info<< "Translating points by " << transVector << endl;

        points += transVector;
    }

    if (args.optionFound("rotate"))
    {
        Pair<vector> n1n2(args.optionLookup("rotate")());
        n1n2[0] /= mag(n1n2[0]);
        n1n2[1] /= mag(n1n2[1]);

        tensor T = rotationTensor(n1n2[0], n1n2[1]);

        Info<< "Rotating points by " << T << endl;

        points = transform(T, points);
    }
    else if (args.optionFound("rollPitchYaw"))
    {
        vector v(args.optionLookup("rollPitchYaw")());

        Info<< "Rotating points by" << nl
            << "    roll  " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    yaw   " << v.z() << endl;


        // Convert to radians
        v *= pi/180.0;

        quaternion R(v.x(), v.y(), v.z());

        Info<< "Rotating points by quaternion " << R << endl;
        points = transform(R, points);
    }
    else if (args.optionFound("yawPitchRoll"))
    {
        vector v(args.optionLookup("yawPitchRoll")());

        Info<< "Rotating points by" << nl
            << "    yaw   " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    roll  " << v.z() << endl;


        // Convert to radians
        v *= pi/180.0;

        scalar yaw = v.x();
        scalar pitch = v.y();
        scalar roll = v.z();

        quaternion R = quaternion(vector(0, 0, 1), yaw);
        R *= quaternion(vector(0, 1, 0), pitch);
        R *= quaternion(vector(1, 0, 0), roll);

        Info<< "Rotating points by quaternion " << R << endl;
        points = transform(R, points);
    }

    if (args.optionFound("scale"))
    {
        vector scaleVector(args.optionLookup("scale")());

        Info<< "Scaling points by " << scaleVector << endl;

        points.replace(vector::X, scaleVector.x()*points.component(vector::X));
        points.replace(vector::Y, scaleVector.y()*points.component(vector::Y));
        points.replace(vector::Z, scaleVector.z()*points.component(vector::Z));
    }

    surf1.movePoints(points);
    surf1.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

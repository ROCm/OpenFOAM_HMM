/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    surfaceTransformPoints

Group
    grpSurfaceUtilities

Description
    Transform (scale/rotate) a surface.
    Like transformPoints but for surfaces.

    The rollPitchYaw and yawPitchRoll options take three angles (degrees)
    that describe the intrinsic Euler rotation.

    rollPitchYaw
    - roll (rotation about X) followed by
    - pitch (rotation about Y) followed by
    - yaw (rotation about Z)

    yawPitchRoll
    - yaw (rotation about Z) followed by
    - pitch (rotation about Y) followed by
    - roll (rotation about X)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Fstream.H"
#include "boundBox.H"
#include "transformField.H"
#include "Pair.H"
#include "Tuple2.H"
#include "axisAngleRotation.H"
#include "EulerCoordinateRotation.H"
#include "MeshedSurfaces.H"

using namespace Foam;
using namespace Foam::coordinateRotations;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transform (translate / rotate / scale) surface points.\n"
        "Like transformPoints but for surfaces.\n"
        "Note: roll=rotate about x, pitch=rotate about y, yaw=rotate about z"
    );
    argList::noParallel();
    argList::addArgument("input", "The input surface file");
    argList::addArgument("output", "The output surface file");
    argList::addOption
    (
        "translate",
        "vector",
        "Translate by specified <vector> - eg, '(1 0 0)' before rotations"
    );
    argList::addOption
    (
        "origin",
        "point",
        "Use specified <point> as origin for rotations"
    );
    argList::addOption
    (
        "rotate",
        "(vectorA vectorB)",
        "Rotate from <vectorA> to <vectorB> - eg, '((1 0 0) (0 0 1))'"
    );
    argList::addOption
    (
        "rotate-angle",
        "(vector scalar)",
        "Rotate <angle> degrees about <vector> - eg, '((1 0 0) 45)'"
    );
    argList::addOption
    (
        "rollPitchYaw",
        "vector",
        "Rotate by '(roll pitch yaw)' degrees"
    );
    argList::addOption
    (
        "yawPitchRoll",
        "vector",
        "Rotate by '(yaw pitch roll)' degrees"
    );
    argList::addOption
    (
        "scale",
        "scalar | vector",
        "Scale by the specified amount - Eg, for uniform [mm] to [m] scaling "
        "use either '(0.001 0.001 0.001)' or simply '0.001'"
    );
    argList args(argc, argv);

    // Verify that an operation has been specified
    {
        const List<word> operationNames
        {
            "translate",
            "rotate",
            "rotate-angle",
            "rollPitchYaw",
            "yawPitchRoll",
            "scale"
        };

        if (!args.count(operationNames))
        {
            FatalError
                << "No operation supplied, "
                << "use at least one of the following:" << nl
                << "   ";

            for (const auto& opName : operationNames)
            {
                FatalError
                    << " -" << opName;
            }

            FatalError
                << nl << exit(FatalError);
        }
    }

    const fileName surfFileName = args[1];
    const fileName outFileName  = args[2];

    Info<< "Reading surf from " << surfFileName << " ..." << nl
        << "Writing surf to " << outFileName << " ..." << endl;

    meshedSurface surf1(surfFileName);

    pointField points(surf1.points());

    vector v;
    if (args.readIfPresent("translate", v))
    {
        Info<< "Translating points by " << v << endl;

        points += v;
    }

    vector origin;
    const bool useOrigin = args.readIfPresent("origin", origin);
    if (useOrigin)
    {
        Info<< "Set origin for rotations to " << origin << endl;
        points -= origin;
    }

    if (args.found("rotate"))
    {
        Pair<vector> n1n2
        (
            args.lookup("rotate")()
        );
        n1n2[0].normalise();
        n1n2[1].normalise();

        const tensor rot(rotationTensor(n1n2[0], n1n2[1]));

        Info<< "Rotating points by " << rot << endl;
        points = transform(rot, points);
    }
    else if (args.found("rotate-angle"))
    {
        const Tuple2<vector, scalar> rotAxisAngle
        (
            args.lookup("rotate-angle")()
        );

        const vector& axis  = rotAxisAngle.first();
        const scalar& angle = rotAxisAngle.second();

        Info<< "Rotating points " << nl
            << "    about " << axis << nl
            << "    angle " << angle << nl;

        const tensor rot(axisAngle::rotation(axis, angle, true));

        Info<< "Rotating points by " << rot << endl;
        points = transform(rot, points);
    }
    else if (args.readIfPresent("rollPitchYaw", v))
    {
        Info<< "Rotating points by" << nl
            << "    roll  " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    yaw   " << v.z() << nl;

        const tensor rot(euler::rotation(euler::eulerOrder::XYZ, v, true));

        Info<< "Rotating points by " << rot << endl;
        points = transform(rot, points);
    }
    else if (args.readIfPresent("yawPitchRoll", v))
    {
        Info<< "Rotating points by" << nl
            << "    yaw   " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    roll  " << v.z() << nl;

        const tensor rot(euler::rotation(euler::eulerOrder::ZYX, v, true));

        Info<< "Rotating points by " << rot << endl;
        points = transform(rot, points);
    }

    List<scalar> scaling;
    if (args.readListIfPresent("scale", scaling))
    {
        // readListIfPresent handles single or multiple values

        if (scaling.size() == 1)
        {
            Info<< "Scaling points uniformly by " << scaling[0] << nl;
            points *= scaling[0];
        }
        else if (scaling.size() == 3)
        {
            Info<< "Scaling points by ("
                << scaling[0] << " "
                << scaling[1] << " "
                << scaling[2] << ")" << nl;

            points.replace(vector::X, scaling[0]*points.component(vector::X));
            points.replace(vector::Y, scaling[1]*points.component(vector::Y));
            points.replace(vector::Z, scaling[2]*points.component(vector::Z));
        }
        else
        {
            FatalError
                << "-scale with 1 or 3 components only" << nl
                << "given: " << args["scale"] << endl
                << exit(FatalError);
        }
    }

    if (useOrigin)
    {
        Info<< "Unset origin for rotations from " << origin << endl;
        points += origin;
    }

    surf1.movePoints(points);
    surf1.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

    The rollPitchYaw option takes three angles (degrees):
    - roll (rotation about x) followed by
    - pitch (rotation about y) followed by
    - yaw (rotation about z)

    The yawPitchRoll does yaw followed by pitch followed by roll.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Fstream.H"
#include "boundBox.H"
#include "transformField.H"
#include "Pair.H"
#include "Tuple2.H"
#include "quaternion.H"
#include "mathematicalConstants.H"

#include "MeshedSurfaces.H"

using namespace Foam;
using namespace Foam::constant::mathematical;


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
    argList::addArgument("surfaceFile");
    argList::addArgument("output surfaceFile");
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

        if (!args.optionCount(operationNames))
        {
            FatalError
                << "No operation supplied, "
                << "use least one of the following:" << nl
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
    if (args.optionReadIfPresent("translate", v))
    {
        Info<< "Translating points by " << v << endl;

        points += v;
    }

    vector origin;
    const bool useOrigin = args.optionReadIfPresent("origin", origin);
    if (useOrigin)
    {
        Info<< "Set origin for rotations to " << origin << endl;
        points -= origin;
    }

    if (args.optionFound("rotate"))
    {
        Pair<vector> n1n2
        (
            args.optionLookup("rotate")()
        );
        n1n2[0] /= mag(n1n2[0]);
        n1n2[1] /= mag(n1n2[1]);

        const tensor rotT = rotationTensor(n1n2[0], n1n2[1]);

        Info<< "Rotating points by " << rotT << endl;

        points = transform(rotT, points);
    }
    else if (args.optionFound("rotate-angle"))
    {
        const Tuple2<vector, scalar> axisAngle
        (
            args.optionLookup("rotate-angle")()
        );

        Info<< "Rotating points " << nl
            << "    about " << axisAngle.first() << nl
            << "    angle " << axisAngle.second() << nl;

        const quaternion quat
        (
            axisAngle.first(),
            axisAngle.second() * pi/180.0  // degToRad
        );

        Info<< "Rotating points by quaternion " << quat << endl;
        points = transform(quat, points);
    }
    else if (args.optionReadIfPresent("rollPitchYaw", v))
    {
        Info<< "Rotating points by" << nl
            << "    roll  " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    yaw   " << v.z() << nl;

        // degToRad
        v *= pi/180.0;

        const quaternion quat(quaternion::rotationSequence::XYZ, v);

        Info<< "Rotating points by quaternion " << quat << endl;
        points = transform(quat, points);
    }
    else if (args.optionReadIfPresent("yawPitchRoll", v))
    {
        Info<< "Rotating points by" << nl
            << "    yaw   " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    roll  " << v.z() << nl;

        // degToRad
        v *= pi/180.0;

        const quaternion quat(quaternion::rotationSequence::ZYX, v);

        Info<< "Rotating points by quaternion " << quat << endl;
        points = transform(quat, points);
    }

    if (args.optionFound("scale"))
    {
        // Use readList to handle single or multiple values
        const List<scalar> scaling = args.optionReadList<scalar>("scale");

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

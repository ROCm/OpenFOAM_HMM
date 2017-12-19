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
    transformPoints

Group
    grpMeshManipulationUtilities

Description
    Transforms the mesh points in the polyMesh directory according to the
    translate, rotate and scale options.

Usage
    Options are:

    -translate vector
        Translates the points by the given vector before rotations

    -rotate (vector vector)
        Rotates the points from the first vector to the second,

    -rotate-angle (vector angle)
        Rotate angle degrees about vector axis.

     or -yawPitchRoll (yawdegrees pitchdegrees rolldegrees)
     or -rollPitchYaw (rolldegrees pitchdegrees yawdegrees)

    -scale scalar|vector
        Scales the points by the given scalar or vector.

    The any or all of the three options may be specified and are processed
    in the above order.

    With -rotateFields (in combination with -rotate/yawPitchRoll/rollPitchYaw)
    it will also read & transform vector & tensor fields.

    Note:
    roll (rotation about x)
    pitch (rotation about y)
    yaw (rotation about z)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ReadFields.H"
#include "pointFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "mathematicalConstants.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void readAndRotateFields
(
    PtrList<GeoField>& flds,
    const fvMesh& mesh,
    const tensor& rotT,
    const IOobjectList& objects
)
{
    ReadFields(mesh, objects, flds);
    forAll(flds, i)
    {
        Info<< "Transforming " << flds[i].name() << endl;
        const dimensionedTensor dimT("t", flds[i].dimensions(), rotT);
        transform(flds[i], dimT, flds[i]);
    }
}


void rotateFields(const argList& args, const Time& runTime, const tensor& T)
{
    #include "createNamedMesh.H"

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    readAndRotateFields(vsFlds, mesh, T, objects);

    PtrList<volVectorField> vvFlds;
    readAndRotateFields(vvFlds, mesh, T, objects);

    PtrList<volSphericalTensorField> vstFlds;
    readAndRotateFields(vstFlds, mesh, T, objects);

    PtrList<volSymmTensorField> vsymtFlds;
    readAndRotateFields(vsymtFlds, mesh, T, objects);

    PtrList<volTensorField> vtFlds;
    readAndRotateFields(vtFlds, mesh, T, objects);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    readAndRotateFields(ssFlds, mesh, T, objects);

    PtrList<surfaceVectorField> svFlds;
    readAndRotateFields(svFlds, mesh, T, objects);

    PtrList<surfaceSphericalTensorField> sstFlds;
    readAndRotateFields(sstFlds, mesh, T, objects);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    readAndRotateFields(ssymtFlds, mesh, T, objects);

    PtrList<surfaceTensorField> stFlds;
    readAndRotateFields(stFlds, mesh, T, objects);

    mesh.write();
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transform (translate / rotate / scale) mesh points.\n"
        "Note: roll=rotate about x, pitch=rotate about y, yaw=rotate about z"
    );
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
    argList::addBoolOption
    (
        "rotateFields",
        "Read and transform vector and tensor fields too"
    );
    argList::addOption
    (
        "scale",
        "scalar | vector",
        "Scale by the specified amount - Eg, for uniform [mm] to [m] scaling "
        "use either '(0.001 0.001 0.001)' or simply '0.001'"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"

    const bool doRotateFields = args.optionFound("rotateFields");

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

    #include "createTime.H"

    word regionName = polyMesh::defaultRegion;
    fileName meshDir = polyMesh::meshSubDir;

    if (args.optionReadIfPresent("region", regionName))
    {
        meshDir = regionName/polyMesh::meshSubDir;
    }

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(meshDir, "points"),
            meshDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

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

        if (doRotateFields)
        {
            rotateFields(args, runTime, rotT);
        }
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

        if (doRotateFields)
        {
            rotateFields(args, runTime, quat.R());
        }
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

        if (doRotateFields)
        {
            rotateFields(args, runTime, quat.R());
        }
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

        if (doRotateFields)
        {
            rotateFields(args, runTime, quat.R());
        }
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


    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< "Writing points into directory " << points.path() << nl << endl;
    points.write();

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //

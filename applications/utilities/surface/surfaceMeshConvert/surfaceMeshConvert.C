/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    surfaceMeshConvert

Group
    grpSurfaceUtilities

Description
    Convert between surface formats with optional scaling or
    transformations (rotate/translate) on a coordinateSystem.

Usage
    \b surfaceMeshConvert inputFile outputFile [OPTION]

    Options:
      - \par -clean
        Perform some surface checking/cleanup on the input surface.

      - \par -scaleIn \<scale\>
        Specify a scaling factor when reading files.

      - \par -scaleOut \<scale\>
        Specify a scaling factor when writing files.

      - \par -dict \<dictionary\>
        Specify an alternative dictionary for constant/coordinateSystems.

      - \par -from \<coordinateSystem\>
        Specify a coordinate System when reading files.

      - \par -to \<coordinateSystem\>
        Specify a coordinate System when writing files.

      - \par -tri
        Triangulate surface.

Note
    The filename extensions are used to determine the file format type.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

#include "MeshedSurfaces.H"
#include "coordinateSystems.H"
#include "cartesianCS.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert between surface formats, using MeshSurface library components"
    );

    argList::noParallel();
    argList::addArgument("input", "The input surface file");
    argList::addArgument("output", "The output surface file");

    argList::addBoolOption
    (
        "clean",
        "Perform some surface checking/cleanup on the input surface"
    );
    argList::addOption
    (
        "scaleIn",
        "factor",
        "Geometry scaling factor on input"
    );
    argList::addOption
    (
        "scaleOut",
        "factor",
        "Geometry scaling factor on output"
    );
    argList::addOption("dict", "file", "Use alternative coordinateSystems");

    argList::addOption
    (
        "from",
        "system",
        "Specify the source coordinate system, applied after '-scaleIn'",
        true // advanced
    );
    argList::addOption
    (
        "to",
        "system",
        "Specify the target coordinate system, applied before '-scaleOut'",
        true // advanced
    );
    argList::addBoolOption
    (
        "tri",
        "Triangulate surface"
    );


    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());

    const fileName importName = args[1];
    const fileName exportName = args[2];

    // disable inplace editing
    if (importName == exportName)
    {
        FatalErrorInFunction
            << "Output file " << exportName << " would overwrite input file."
            << exit(FatalError);
    }

    // Check that reading/writing is supported
    if
    (
        !MeshedSurface<face>::canRead(importName, true)
     || !MeshedSurface<face>::canWriteType(exportName.ext(), true)
    )
    {
        return 1;
    }


    // The coordinate transformations (must be cartesian)
    autoPtr<coordSystem::cartesian> fromCsys;
    autoPtr<coordSystem::cartesian> toCsys;

    if (args.found("from") || args.found("to"))
    {
        IOobject ioCsys = IOobject::selectIO
        (
            IOobject
            (
                coordinateSystems::typeName,
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            args.get<fileName>("dict", "")
        );

        if (!ioCsys.typeHeaderOk<coordinateSystems>(false))
        {
            FatalErrorInFunction
                << "Cannot open coordinateSystems file\n    "
                << ioCsys.objectPath() << nl
                << exit(FatalError);
        }

        coordinateSystems globalCoords(ioCsys);

        if (args.found("from"))
        {
            const word csName(args["from"]);
            const auto* csPtr = globalCoords.lookupPtr(csName);

            if (!csPtr)
            {
                FatalErrorInFunction
                    << "Cannot find -from " << csName << nl
                    << "available coordinateSystems: "
                    << flatOutput(globalCoords.names()) << nl
                    << exit(FatalError);
            }

            fromCsys = autoPtr<coordSystem::cartesian>::New(*csPtr);
        }

        if (args.found("to"))
        {
            const word csName(args["to"]);
            const auto* csPtr = globalCoords.lookupPtr(csName);

            if (!csPtr)
            {
                FatalErrorInFunction
                    << "Cannot find -to " << csName << nl
                    << "available coordinateSystems: "
                    << flatOutput(globalCoords.names()) << nl
                    << exit(FatalError);
            }

            toCsys = autoPtr<coordSystem::cartesian>::New(*csPtr);
        }

        // Maybe fix this later
        if (fromCsys && toCsys)
        {
            FatalErrorInFunction
                << "Only allowed '-from' or '-to' option at the moment."
                << exit(FatalError);
        }
    }


    {
        MeshedSurface<face> surf(importName);

        if (args.found("clean"))
        {
            surf.cleanup(true);
        }

        scalar scaleIn = 0;
        if (args.readIfPresent("scaleIn", scaleIn) && scaleIn > 0)
        {
            Info<< "scale input " << scaleIn << endl;
            surf.scalePoints(scaleIn);
        }

        if (fromCsys)
        {
            Info<< "move points from coordinate system: "
                << fromCsys->name() << endl;
            tmp<pointField> tpf = fromCsys->localPosition(surf.points());
            surf.movePoints(tpf());
        }

        if (toCsys)
        {
            Info<< "move points to coordinate system: "
                << toCsys->name() << endl;
            tmp<pointField> tpf = toCsys->globalPosition(surf.points());
            surf.movePoints(tpf());
        }

        scalar scaleOut = 0;
        if (args.readIfPresent("scaleOut", scaleOut) && scaleOut > 0)
        {
            Info<< "scale output " << scaleOut << endl;
            surf.scalePoints(scaleOut);
        }

        if (args.found("tri"))
        {
            Info<< "triangulate" << endl;
            surf.triangulate();
        }

        Info<< "writing " << exportName;
        surf.write(exportName);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //

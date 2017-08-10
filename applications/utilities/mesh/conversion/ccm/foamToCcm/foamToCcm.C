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
    foamToCcm

Group
    grpMeshConversionUtilities

Description
    Translates OPENFOAM mesh and/or results to CCM format

Usage
    \b foamToCcm [OPTION]

    Options:
      - \par -mesh
        convert mesh only to CCM format

      - \par -name \<name\>
        Provide alternative base name. Default is <tt>meshExport</tt>.

      - \par -overwrite
        No backup of existing output files.

      - \par -remap \<name\>
        Use specified remapping dictionary instead of
        <tt>constant/remapping</tt>

      - \par -results
        Convert results only to CCM format

Note
    - No parallel data
    - No Lagrangian elements
    - the -noZero time option can be useful to avoid the often incomplete
      initial conditions (missing useful calculated values)

See also
    Foam::ccm::writer for information about the
    <tt>constant/remapping</tt> file.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "volFields.H"
#include "OFstream.H"
#include "IOobjectList.H"
#include "scalarIOField.H"
#include "tensorIOField.H"

#include "ccm.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Translate OPENFOAM data to CCM format"
    );

    Foam::timeSelector::addOptions();
    argList::noParallel();
    argList::addBoolOption
    (
        "mesh",
        "convert mesh only"
    );
    argList::addOption
    (
        "name",
        "name",
        "provide alternative base name. Default is <meshExport>."
    );
    argList::addBoolOption
    (
        "overwrite",
        "no backup of existing output files"
    );
    argList::addOption
    (
        "remap",
        "name",
        "use specified remapping dictionary instead of <constant/remapping>"
    );
    argList::addBoolOption
    (
        "results",
        "convert results only"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    // get times list
    instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    const bool optMesh      = args.optionFound("mesh");
    const bool optResults   = args.optionFound("results");
    const bool optOverwrite = args.optionFound("overwrite");

    fileName exportName = ccm::writer::defaultMeshName;
    if (args.optionReadIfPresent("name", exportName))
    {
        const word ext = exportName.ext();
        // strip erroneous extension (.ccm, .ccmg, .ccmp)
        if (ext == "ccm" || ext == "ccmg" || ext == "ccmp")
        {
            exportName = exportName.lessExt();
        }
    }
    else if (args.optionFound("case"))
    {
        exportName += '-' + args.globalCaseName();
    }

    if (optMesh && optResults)
    {
        Warning
            << "\n-mesh and -results options are mutually exclusive\n"
            << endl;
        args.printUsage();
        FatalError.exit();
    }

//     // skip over time=0, unless some other time option has been specified
//     if
//     (
//         !args.optionFound("zeroTime")
//      && !args.optionFound("time")
//      && !args.optionFound("latestTime")
//      && Times.size() > 2
//     )
//     {
//         startTime = 2;
//     }
//
//    runTime.setTime(Times[startTime], startTime);

    runTime.setTime(timeDirs[0], 0);
    if (optMesh)
    {
        // convert mesh only
        #include "createPolyMesh.H"

        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);

            #include "getTimeIndex.H"

            if (timeI == 0)
            {
                ccm::writer writer
                (
                    exportName + ".ccmg",
                    mesh,
                    !optOverwrite
                );
                writer.writeGeometry();
            }
            else if (mesh.moving())
            {
                ccm::writer writer
                (
                    exportName + ".ccmg_" + timeName,
                    mesh,
                    !optOverwrite
                );
                writer.writeGeometry();
            }
        }
    }
    else
    {
        // convert fields with or without converting mesh
        #include "createMesh.H"

        // #include "checkHasMovingMesh.H"
        // #include "checkHasLagrangian.H"

        IOobjectList objects(mesh, timeDirs.last().name());

        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);

            #include "getTimeIndex.H"

            Info<< "has "
                << mesh.nCells() << " cells, "
                << mesh.nPoints() << " points, "
                << mesh.boundaryMesh().size() << " patches"
                << endl;

            if (!optResults)
            {
                if (timeI == 0)
                {
                    ccm::writer writer
                    (
                        exportName + ".ccmg",
                        mesh,
                        !optOverwrite
                    );
                    writer.writeGeometry();
                }
                else if (mesh.moving())
                {
                    ccm::writer writer
                    (
                        exportName + ".ccmg_" + timeName,
                        mesh,
                        !optOverwrite
                    );
                    writer.writeGeometry();
                }
            }

            ccm::writer writer
            (
                exportName + ".ccmp_" + timeName,
                mesh,
                !optOverwrite
            );
            // writer.setTopologyFile(exportName + ".ccmg");
            Info<< "writing solution:";
            if (args.optionFound("remap"))
            {
                writer.writeSolution(objects, args["remap"]);
            }
            else
            {
                writer.writeSolution(objects);
            }
        }
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //

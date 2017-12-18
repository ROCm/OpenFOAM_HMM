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
    ccmToFoam

Group
    grpMeshConversionUtilities

Description
    Reads CCM files as written by PROSTAR/STARCCM and writes an
    OPENFOAM polyMesh.

Usage
    \b ccmToFoam [OPTION] ccmMesh

    Options:
      - \par -ascii
        Write in ASCII format instead of binary

      - \par -export
        re-export mesh in CCM format for post-processing

      - \par -list
        List some information about the geometry

      - \par -name \<name\>
        Provide alternative base name for export.
        Default is <tt>meshExport</tt>.

      - \par -noBaffles
        Remove any baffles by merging the faces.

      - \par -merge
        Merge in-place interfaces

      - \par -numbered
        Use numbered patch/zone (not names) directly from ccm ids.

      - \par -remap \<name\>
        Use specified remapping dictionary instead of
        <tt>constant/remapping</tt>

      - \par -scale \<factor\>
        Specify an alternative geometry scaling factor.
        The default is \b 1 (no scaling).

      - \par -solids
        Treat any solid cells present just like fluid cells.
        The default is to remove them.

Note
    - sub-domains (fluid | solid | porosity) are stored as separate domains
      within the CCM file. These are merged together to form a single mesh.
    - baffles are written as interfaces for later use

See also
    Foam::ccm::reader for more information about the File Locations

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "ccm.H"
#include "regionSplit.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Reads CCM files as written by PROSTAR/STARCCM and writes an OPENFOAM "
        " polyMesh. Multi-region support for PROSTAR meshes should be stable."
        " Multi-region merging for STARCCM meshes will not always be"
        " successful."
    );

    argList::noParallel();
    argList::addArgument("ccmMesh");
    argList::addBoolOption
    (
        "ascii",
        "write in ASCII format instead of binary"
    );
    argList::addBoolOption
    (
        "export",
        "re-export mesh in CCM format for post-processing"
    );
    argList::addBoolOption
    (
        "list",
        "list some information about the geometry"
    );
    argList::addOption
    (
        "remap",
        "name",
        "use specified remapping dictionary instead of <constant/remapping>"
    );
    argList::addOption
    (
        "name",
        "name",
        "provide alternative base name when re-exporting (implies -export). "
        "Default is <meshExport>."
    );
    argList::addBoolOption
    (
        "noBaffles",
        "remove any baffles by merging the faces"
    );
    argList::addBoolOption
    (
        "merge",
        "merge in-place interfaces"
    );
    argList::addBoolOption
    (
        "numbered",
        "use numbered names (eg, patch_0, zone_0) only"
    );
    argList::addOption
    (
        "scale",
        "scale",
        "geometry scaling factor - default is 1 (ie, no scaling)"
    );
    argList::addBoolOption
    (
        "solids",
        "treat any solid cells present just like fluid cells. "
        "the default is to remove them."
    );

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());
    runTime.functionObjects().off();

    const bool optList = args.optionFound("list");

    // exportName only has a size when export is in effect
    fileName exportName;
    if (args.optionReadIfPresent("name", exportName))
    {
        const word ext = exportName.ext();
        // strip erroneous extension (.ccm, .ccmg, .ccmp)
        if (ext == "ccm" || ext == "ccmg" || ext == "ccmp")
        {
            exportName = exportName.lessExt();
        }
    }
    else if (args.optionFound("export"))
    {
        exportName = ccm::writer::defaultMeshName;
        if (args.optionFound("case"))
        {
            exportName += '-' + args.globalCaseName();
        }
    }

    // By default, no scaling
    const scalar scaleFactor = args.optionLookupOrDefault("scale", 1.0);

    // Default to binary output, unless otherwise specified
    const IOstream::streamFormat format =
    (
        args.optionFound("ascii")
      ? IOstream::ASCII
      : IOstream::BINARY
    );

    // Increase the precision of the points data
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));


    // Read control options
    // ~~~~~~~~~~~~~~~~~~~~

    ccm::reader::options rOpts;
    rOpts.removeBaffles(args.optionFound("noBaffles"));
    rOpts.mergeInterfaces(args.optionFound("merge"));

    if (args.optionFound("numbered"))
    {
        rOpts.useNumberedNames(true);
    }

    if (args.optionFound("solids"))
    {
        Info<< "treating solids like fluids" << endl;
        rOpts.keepSolid(true);
    }
    else
    {
        rOpts.keepSolid(false);
    }

    // CCM reader for reading geometry/solution
    ccm::reader reader(args[1], rOpts);

    // list the geometry information
    if (optList)
    {
        Info<< "mesh geometry information:" << endl;
        if (reader.hasGeometry())
        {
            Info<< nl << "cellTable:" << reader.cellTableInfo()
                << nl << "boundaryRegion:" << reader.boundaryTableInfo()
                << nl << "interfaces:" << reader.interfaceDefinitionsInfo()
                << endl;

            if
            (
                args.optionFound("remap")
              ? reader.remapMeshInfo(runTime, args["remap"])
              : reader.remapMeshInfo(runTime)
            )
            {
                Info<< nl
                    << "Remapped cellTable:" << reader.cellTableInfo() << nl
                    << "Remapped boundaryRegion:" << reader.boundaryTableInfo()
                    << endl;
            }
        }
        else
        {
            Info<< "NONE" << endl;
        }

        return 0;
    }
    else if (reader.readGeometry(scaleFactor))
    {
        autoPtr<polyMesh> mesh =
        (
            args.optionFound("remap")
          ? reader.mesh(runTime, args["remap"])
          : reader.mesh(runTime)
        );

        // report mesh bounding box information
        Info<< nl << "Bounding box size: " << mesh().bounds().span() << nl;

        // check number of regions
        regionSplit rs(mesh);

        Info<< "Number of regions: " << rs.nRegions();
        if (rs.nRegions() == 1)
        {
            Info<< " (OK)." << nl;
        }
        else
        {
            Info<< nl << nl
                << "**************************************************" << nl
                << "**  WARNING: the mesh has disconnected regions  **" << nl
                << "**************************************************" << nl;
        }
        Info<< endl;
        reader.writeMesh(mesh, format);

        // exportName only has a size when export is in effect
        if (exportName.size())
        {
            const fileName geomName = exportName + ".ccmg";
            Info<< nl << "Re-exporting geometry as " << geomName << nl;
            ccm::writer(geomName, mesh()).writeGeometry();
        }
    }
    else
    {
        FatalErrorInFunction
            << "could not read geometry"
            << exit(FatalError);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Code chunk for post-processing conversion of cloud(s) to VTK PolyData
    (.vtp extension).

\*---------------------------------------------------------------------------*/

if (doLagrangian)
{
    const fileName cloudPrefix = (regionDir/cloud::prefix);

    wordList cloudNames = ListOps::create<word>
    (
        readDir
        (
            mesh.time().path()/mesh.time().timeName()/cloudPrefix,
            fileName::DIRECTORY
        ),
        nameOp<fileName>()
    );

    // Synchronise cloud names
    Pstream::combineGather(cloudNames, ListOps::uniqueEqOp<word>());
    Pstream::combineScatter(cloudNames);

    // Consistent order
    Foam::sort(cloudNames);


    for (const word& cloudName : cloudNames)
    {
        IOobjectList cloudObjs(mesh, runTime.timeName(), cloudPrefix/cloudName);

        bool isCloud = false;
        if (cloudObjs.erase("coordinates"))
        {
            isCloud = true;
        }
        if (cloudObjs.erase("positions"))
        {
            isCloud = true;
        }

        if (!returnReduce(isCloud, orOp<bool>()))
        {
            continue;
        }

        // Limited to basic IOField types
        cloudObjs.filterClasses
        (
            stringListOps::foundOp<word>(fieldTypes::basic)
        );

        // Are there cloud fields (globally)?
        if (returnReduce(cloudObjs.empty(), andOp<bool>()))
        {
            continue;
        }

        vtk::lagrangianWriter writer
        (
            meshProxy.baseMesh(),
            cloudName,
            writeOpts,
            // Output name for the cloud
            (
                outputDir/regionDir/cloud::prefix
              / cloudName/cloudName + timeDesc
            ),
            Pstream::parRun()
        );

        Info<< "    Lagrangian: "
            << args.relativePath(writer.output()) << nl;

        writer.writeTimeValue(mesh.time().value());
        writer.writeGeometry();

        // Begin CellData/PointData
        writer.beginParcelData();

        writer.writeFields<label>(cloudObjs);
        writer.writeFields<scalar>(cloudObjs);
        writer.writeFields<vector>(cloudObjs);
        writer.writeFields<sphericalTensor>(cloudObjs);
        writer.writeFields<symmTensor>(cloudObjs);
        writer.writeFields<tensor>(cloudObjs);

        fileName outputName(writer.output());

        writer.close();

        if (Pstream::master())
        {
            // Add to file-series and emit as JSON

            fileName seriesName(vtk::seriesWriter::base(outputName));

            vtk::seriesWriter& series = vtkSeries(seriesName);

            // First time?
            // Load from file, verify against filesystem,
            // prune time >= currentTime
            if (series.empty())
            {
                series.load(seriesName, true, timeValue);
            }

            series.append(timeValue, outputName);
            series.write(seriesName);
        }
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Code chunk for converting finite-area - included by foamToVTK.

    Typedef vtkWriterType_areaMesh is declared in writeAreaFields.H

\*---------------------------------------------------------------------------*/

//
// Finite-area mesh and fields - need not exist
//

// No subsetting!
if (doFiniteArea)
{
    using reportFields = foamToVtkReportFields;

    autoPtr<faMesh> faMeshPtr;

    const label nAreaFields =
        objects.count(stringListOps::foundOp<word>(fieldTypes::area));

    if (nAreaFields || withMeshIds)
    {
        faMeshPtr = faMesh::TryNew(meshProxy.baseMesh());
    }

    if (faMeshPtr && (nAreaFields || withMeshIds))
    {
        const faMesh& areaMesh = faMeshPtr();

        reportFields::area(Info, objects);

        const auto& pp = faMeshPtr->patch();

        vtkWriterType_areaMesh writer
        (
            pp,
            writeOpts,
            (
                outputDir/regionDir/"finite-area"
              / "finiteArea" + timeDesc
            ),
            Pstream::parRun()
        );
        Info<< "    Area      : "
            << args.relativePath(writer.output()) << nl;

        writer.beginFile(areaMesh.name());

        writer.writeTimeValue(timeValue);
        writer.writeGeometry();

        // Optionally with (cellID, faceLabels, procID) fields
        writer.beginCellData
        (
            (withMeshIds ? 2 + (writer.parallel() ? 1 : 0) : 0)
          + nAreaFields
        );

        if (withMeshIds)
        {
            const globalIndex procAddr(areaMesh.nFaces());

            // Use global indexed values for the 'cell' ids
            writer.writeCellData("cellID", identity(procAddr.range()));

            // Use proc-local data for faceLabels
            // (confusing enough already without renumbering)
            writer.writeCellData("faceLabels", areaMesh.faceLabels());

            writer.writeProcIDs();  // parallel only
        }

        writeAllAreaFields
        (
            writer,
            areaMesh,
            objects,
            true // syncPar
        );

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

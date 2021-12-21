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

    if (nAreaFields)
    {
        const bool oldThrowingError = FatalError.throwing(true);
        try
        {
            faMeshPtr.reset(new faMesh(meshProxy.baseMesh()));
        }
        catch (const Foam::error& err)
        {
            faMeshPtr.clear();
        }
        FatalError.throwing(oldThrowingError);
    }

    if (faMeshPtr && nAreaFields)
    {
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

        writer.beginFile(faMeshPtr->name());

        writer.writeTimeValue(timeValue);
        writer.writeGeometry();

        writer.beginCellData(nAreaFields);

        writeAllAreaFields
        (
            writer,
            *faMeshPtr,
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

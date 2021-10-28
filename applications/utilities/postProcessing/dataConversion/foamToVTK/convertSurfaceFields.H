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
    Code chunk for post-processing surface fields to VTK PolyData

\*---------------------------------------------------------------------------*/

{
    using reportFields = foamToVtkReportFields;

    // Load only once when possible
    label nSurfaceScalarField = -1;
    label nSurfaceVectorField = -1;

    PtrList<const surfaceScalarField> sScalars;
    PtrList<const surfaceVectorField> sVectors;

    // Surface Fields
    if (doSurfaceFields)
    {
        if (nSurfaceScalarField == -1)
        {
            sScalars = readFields<surfaceScalarField>
            (
                meshProxy,
                objects,
                selectedFields
            );

            reportFields::print("    surfScalar   :", Info, sScalars);
            nSurfaceScalarField = sScalars.size();
        }
        else
        {
            sScalars.resize(nSurfaceScalarField); // Consistent sizing
        }

        if (nSurfaceVectorField == -1)
        {
            sVectors = readFields<surfaceVectorField>
            (
                meshProxy,
                objects,
                selectedFields
            );

            reportFields::print("    surfVector   :", Info, sVectors);
            nSurfaceVectorField = sVectors.size();
        }
        else
        {
            sVectors.resize(nSurfaceVectorField); // Consistent sizing
        }

        if (sScalars.size())
        {
            // Change scalar fields into vector fields, but leave
            // the count of vector fields unchanged. This allows us to
            // easily delete these synthetic fields later.

            surfaceVectorField unitNorm(mesh.Sf()/mesh.magSf());

            sVectors.resize(nSurfaceVectorField + nSurfaceScalarField);

            label nExtra = 0;
            for (const auto& ssf : sScalars)
            {
                surfaceVectorField* tsvfPtr = (ssf * unitNorm).ptr();
                tsvfPtr->rename(ssf.name());
                sVectors.set(nSurfaceVectorField + nExtra, tsvfPtr);
                ++nExtra;
            }
        }

        if (sVectors.size())
        {
            vtk::surfaceFieldWriter writer
            (
                meshProxy.mesh(),
                writeOpts,
                (
                    outputDir/regionDir
                  / "surface-fields"/"surfaceFields" + timeDesc
                ),
                Pstream::parRun()
            );

            Info<< "    Surface   : "
                << args.relativePath(writer.output()) << nl;


            writer.writeTimeValue(timeValue);
            writer.writeGeometry();

            writer.beginPointData(sVectors.size());

            for (const auto& fld : sVectors)
            {
                writer.write(fld);
            }

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


    // Write faceZones (POLYDATA file, one for each zone)

    if (!selectedFaceZones.empty() && !mesh.faceZones().empty())
    {
        if (nSurfaceScalarField == -1)
        {
            sScalars = readFields<surfaceScalarField>
            (
                meshProxy,
                objects,
                selectedFields
            );
            nSurfaceScalarField = sScalars.size();

            reportFields::print("    surfScalar   :", Info, sScalars);
        }
        else
        {
            sScalars.resize(nSurfaceScalarField); // Consistent sizing
        }

        if (nSurfaceVectorField == -1)
        {
            sVectors = readFields<surfaceVectorField>
            (
                meshProxy,
                objects,
                selectedFields
            );
            nSurfaceVectorField = sVectors.size();

            reportFields::print("    surfVector   :", Info, sVectors);
        }
        else
        {
            sVectors.resize(nSurfaceVectorField); // Consistent sizing
        }

        for (const faceZone& fz : mesh.faceZones())
        {
            if (!selectedFaceZones.match(fz.name()))
            {
                continue;
            }

            // Retrieve as primitiveFacePatch with faces properly flipped
            const primitiveFacePatch& pp = fz();

            vtkWriterType_faceZone writer
            (
                pp,
                writeOpts,
                (
                    outputDir/regionDir/fz.name()
                  / (meshProxy.useSubMesh() ? meshProxy.name() : fz.name())
                  + timeDesc
                ),
                Pstream::parRun()
            );

            Info<< "    FaceZone  : "
                << args.relativePath(writer.output()) << nl;


            writer.beginFile(fz.name());
            writer.writeTimeValue(timeValue);
            writer.writeGeometry();

            writer.beginCellData(sScalars.size() + sVectors.size());

            for (const auto& fld : sScalars)
            {
                writer.write(fld, fz.addressing());
            }
            for (const auto& fld : sVectors)
            {
                writer.write(fld, fz.addressing());
            }

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
}


// ************************************************************************* //

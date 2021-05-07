/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Code chunk for converting volume fields on processor boundaries,
    included by foamToVTK.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Generate processor/processorN/procBoundary

{
    using reportFields = foamToVtkReportFields;

    const label nVolFields =
    (
        objects.count(stringListOps::foundOp<word>(fieldTypes::volume))
    );

    reportFields::volume(Info, objects);

    // Setup for the vtm writer.

    fileName vtmOutputBase
    (
        outputDir/regionDir/vtkName + timeDesc
    );

    // Naming
    const auto subDirNaming =
        [](const label i) -> fileName
        { return "processor" / ("processor" + Foam::name(i)); };


    // Dummy writer.
    autoPtr<vtk::internalWriter> internalWriter;

    // Setup the patch writers
    PtrList<vtk::patchWriter> patchWriters;

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelList patchIds =
        identity
        (
            patches.size()-patches.nNonProcessor(),
            patches.nNonProcessor()
        );

    forAll(patchIds, i)
    {
        if (!isA<processorPolyPatch>(patches[patchIds[i]]))
        {
            patchIds.resize(i);
            break;
        }
    }

    patchWriters.resize(patchIds.size());

    label nPatchWriters = 0;

    List<wordList> procPatchNames(Pstream::nProcs());

    procPatchNames[Pstream::myProcNo()].resize(patchIds.size());

    for (const label patchId : patchIds)
    {
        const polyPatch& pp = patches[patchId];

        auto writer = autoPtr<vtk::patchWriter>::New
        (
            meshProxy.mesh(),
            labelList(one{}, pp.index()),
            writeOpts,
            nearCellValue,
            (
                vtmOutputBase
              / subDirNaming(Pstream::myProcNo())
              / pp.name()
            ),
            false  // This MUST be non-parallel (serial only)
        );

        procPatchNames[Pstream::myProcNo()][nPatchWriters] = pp.name();

        writer->writeTimeValue(timeValue);
        writer->writeGeometry();

        // Transfer writer to list for later use
        patchWriters.set(nPatchWriters++, writer);
    }
    patchWriters.resize(nPatchWriters);


    Pstream::gatherList(procPatchNames);

    // CellData
    {
        for (vtk::patchWriter& writer : patchWriters)
        {
            // Optionally with patchID, procID, neighID fields
            // - use Pstream::parRun() not writer.parallel() !!
            writer.beginCellData
            (
                (withMeshIds ? 1 + (Pstream::parRun() ? 2 : 0) : 0)
              + nVolFields
            );

            if (withMeshIds)
            {
                writer.writePatchIDs();
                writer.writeProcIDs();  // parallel only
                writer.writeNeighIDs(); // parallel only
            }
        }

        writeAllVolFields
        (
            internalWriter,
            patchWriters,
            meshProxy,
            objects,
            true  // syncPar
        );

        // End CellData is implicit
    }


    // Finish writers
    if (internalWriter)
    {
        internalWriter->close();
    }

    for (vtk::patchWriter& writer : patchWriters)
    {
        writer.close();
    }

    patchWriters.clear();


    // Collective output

    const label nProcPatches = returnReduce(nPatchWriters, sumOp<label>());

    if (Pstream::master() && nProcPatches)
    {
        Info<< "Wrote " << nProcPatches << " processor boundaries from "
            << Pstream::nProcs() << " processes" << nl;


        // Collect individual boundaries into a vtm file
        vtk::vtmWriter vtmBoundaries;

        // Naming for vtm
        fileName outputName(vtmOutputBase / "processor");
        outputName.ext(vtmBoundaries.ext());

        vtmBoundaries.setTime(timeValue);

        forAll(procPatchNames, proci)
        {
            label n = 0;

            const word blockName("proc" + Foam::name(proci));
            const fileName dirName = subDirNaming(proci);

            for (const word& patchName : procPatchNames[proci])
            {
                if (!n)
                {
                    vtmBoundaries.beginBlock(blockName);
                    ++n;
                }

                vtmBoundaries.append_vtp
                (
                    patchName,
                    dirName/patchName
                );
            }

            if (n)
            {
                vtmBoundaries.endBlock();
            }
        }


        // Emit "processor.vtm" with collection of processor boundaries
        vtmBoundaries.write(outputName);
    }
}


// ************************************************************************* //

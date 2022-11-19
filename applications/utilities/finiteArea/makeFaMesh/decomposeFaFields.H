/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Write proc addressing and decompose area fields (parallel only).

\*---------------------------------------------------------------------------*/

// Embed do-while to support early termination
if (doDecompose && Pstream::parRun())
do
{
    faMeshReconstructor reconstructor(aMesh, IOobjectOption::READ_IF_PRESENT);

    if (!reconstructor.good())
    {
        Info<< "Missing volume proc-addressing, "
               "cannot generate area proc-addressing." << nl
            << "Also skip decomposing area fields...."
            << endl;
        break;
    }

    reconstructor.writeAddressing();

    Info<< "Wrote proc-addressing" << nl << endl;

    // Handle area fields
    // ------------------

    faFieldDecomposer::fieldsCache areaFieldsCache;

    const faMesh& fullMesh = reconstructor.mesh();

    {
        // Use uncollated (or master uncollated) file handler here.
        // - each processor is reading in the identical serial fields.
        // - nothing should be parallel-coordinated.

        // Similarly, if we write the serial finite-area mesh, this is only
        // done from one processor!

        reconstructor.writeMesh();

        if (doDecompFields)
        {
            const bool oldDistributed = fileHandler().distributed();
            auto oldHandler = fileHandler(fileOperation::NewUncollated());
            fileHandler().distributed(true);

            IOobjectList objects(fullMesh.time(), runTime.timeName());

            areaFieldsCache.readAllFields(fullMesh, objects);

            // Restore old settings
            if (oldHandler)
            {
                fileHandler(std::move(oldHandler));
            }
            fileHandler().distributed(oldDistributed);
        }
    }

    const label nAreaFields = areaFieldsCache.size();

    if (nAreaFields)
    {
        Info<< "Decomposing " << nAreaFields << " area fields" << nl;

        faFieldDecomposer fieldDecomposer
        (
            fullMesh,
            aMesh,
            reconstructor.edgeProcAddressing(),
            reconstructor.faceProcAddressing(),
            reconstructor.boundaryProcAddressing()
        );

        // Report
        areaFieldsCache.decomposeAllFields(fieldDecomposer, true);
        Info<< endl;
    }
} while (false);


// ************************************************************************* //

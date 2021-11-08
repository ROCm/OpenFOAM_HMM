/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Write proc addressing and decompose area fields (parallel only).

\*---------------------------------------------------------------------------*/

if (doDecompose && Pstream::parRun())
{
    faMeshReconstructor reconstructor(aMesh);
    reconstructor.writeAddressing();

    Info<< "Wrote proc-addressing" << nl << endl;

    // Handle area fields
    // ------------------

    PtrList<areaScalarField> areaScalarFields;
    PtrList<areaVectorField> areaVectorFields;
    PtrList<areaSphericalTensorField> areaSphTensorFields;
    PtrList<areaSymmTensorField> areaSymmTensorFields;
    PtrList<areaTensorField> areaTensorFields;

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

            faFieldDecomposer::readFields
            (
                fullMesh, objects, areaScalarFields
            );
            faFieldDecomposer::readFields
            (
                fullMesh, objects, areaVectorFields
            );
            faFieldDecomposer::readFields
            (
                fullMesh, objects, areaSphTensorFields
            );
            faFieldDecomposer::readFields
            (
                fullMesh, objects, areaSymmTensorFields
            );
            faFieldDecomposer::readFields
            (
                fullMesh, objects, areaTensorFields
            );

            // Restore old settings
            if (oldHandler)
            {
                fileHandler(std::move(oldHandler));
            }
            fileHandler().distributed(oldDistributed);
        }
    }

    const label nAreaFields =
    (
        areaScalarFields.size()
      + areaVectorFields.size()
      + areaSphTensorFields.size()
      + areaSymmTensorFields.size()
      + areaTensorFields.size()
    );

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

        if (areaScalarFields.size())
        {
            Info<< "  scalars: "
                << flatOutput(PtrListOps::names(areaScalarFields)) << nl;
            fieldDecomposer.decomposeFields(areaScalarFields);
        }

        if (areaVectorFields.size())
        {
            Info<< "  vectors: "
                << flatOutput(PtrListOps::names(areaVectorFields)) << nl;
            fieldDecomposer.decomposeFields(areaVectorFields);
        }

        if (areaSphTensorFields.size())
        {
            Info<< "  sphTensors: "
                << flatOutput(PtrListOps::names(areaSphTensorFields)) << nl;
            fieldDecomposer.decomposeFields(areaSphTensorFields);
        }

        if (areaSymmTensorFields.size())
        {
            Info<< "  symmTensors: "
                << flatOutput(PtrListOps::names(areaSymmTensorFields)) << nl;
            fieldDecomposer.decomposeFields(areaSymmTensorFields);
        }

        if (areaTensorFields.size())
        {
            Info<< "  tensors: "
                << flatOutput(PtrListOps::names(areaTensorFields)) << nl;
            fieldDecomposer.decomposeFields(areaTensorFields);
        }

        Info<< endl;
    }
}


// ************************************************************************* //

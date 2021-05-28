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
    Decompose area fields, when mesh was generated in parallel

\*---------------------------------------------------------------------------*/

if (Pstream::parRun())
{
    faMeshReconstructor reconstructor(areaMesh);
    reconstructor.writeAddressing();

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

        const bool oldDistributed = fileHandler().distributed();
        auto oldHandler = fileHandler(fileOperation::NewUncollated());
        fileHandler().distributed(true);

        IOobjectList objects(fullMesh.time(), runTime.timeName());

        faFieldDecomposer::readFields(fullMesh, objects, areaScalarFields);
        faFieldDecomposer::readFields(fullMesh, objects, areaVectorFields);
        faFieldDecomposer::readFields(fullMesh, objects, areaSphTensorFields);
        faFieldDecomposer::readFields(fullMesh, objects, areaSymmTensorFields);
        faFieldDecomposer::readFields(fullMesh, objects, areaTensorFields);

        // Restore old settings
        if (oldHandler)
        {
            fileHandler(std::move(oldHandler));
        }
        fileHandler().distributed(oldDistributed);
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
        Info<< "Decomposing " << nAreaFields << " area fields" << nl << endl;

        faFieldDecomposer fieldDecomposer
        (
            fullMesh,
            areaMesh,
            reconstructor.edgeProcAddressing(),
            reconstructor.faceProcAddressing(),
            reconstructor.boundaryProcAddressing()
        );

        fieldDecomposer.decomposeFields(areaScalarFields);
        fieldDecomposer.decomposeFields(areaVectorFields);
        fieldDecomposer.decomposeFields(areaSphTensorFields);
        fieldDecomposer.decomposeFields(areaSymmTensorFields);
        fieldDecomposer.decomposeFields(areaTensorFields);
    }
}

// ************************************************************************* //

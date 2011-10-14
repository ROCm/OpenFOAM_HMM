//
// createAllMesh.H
// ~~~~~~~~~~~~~~~

    autoPtr<fvMesh> allMesh;
    autoPtr<fvFieldReconstructor> procToAllMapper;
    PtrList<fvFieldDecomposer> allToProcMappers;

    if (temperatureCoupled)
    {
        Foam::Info
            << "Create mesh for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;

        allMesh.reset
        (
            new Foam::fvMesh
            (
                Foam::IOobject
                (
                    Foam::fvMesh::defaultRegion,
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                )
            )
        );

        procToAllMapper.reset
        (
            new fvFieldReconstructor
            (
                allMesh(),
                procMeshes,
                faceProcAddressing,
                cellProcAddressing,
                boundaryProcAddressing
            )
        );


        allToProcMappers.setSize
        (
            rp.fluidRegionNames().size()
          + rp.solidRegionNames().size()
        );

        forAll(allToProcMappers, i)
        {
            allToProcMappers.set
            (
                i,
                new fvFieldDecomposer
                (
                    allMesh(),
                    procMeshes[i],
                    faceProcAddressing[i],
                    cellProcAddressing[i],
                    boundaryProcAddressing[i]
                )
            );
        }
    }

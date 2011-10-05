    forAll(rp.fluidRegionNames(), i)
    {
        Info<< "Create fluid mesh for region " << rp.fluidRegionNames()[i]
            << " for time = " << runTime.timeName() << nl << endl;

        label procI = fluidToProc[i];

        procMeshes.set
        (
            procI,
            new fvMesh
            (
                IOobject
                (
                    rp.fluidRegionNames()[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

        if (temperatureCoupled)
        {
            cellProcAddressing.set
            (
                procI,
                new labelIOList
                (
                    IOobject
                    (
                        "cellRegionAddressing",
                        procMeshes[procI].facesInstance(),
                        procMeshes[procI].meshSubDir,
                        procMeshes[procI],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );

            faceProcAddressing.set
            (
                procI,
                new labelIOList
                (
                    IOobject
                    (
                        "faceRegionAddressing",
                        procMeshes[procI].facesInstance(),
                        procMeshes[procI].meshSubDir,
                        procMeshes[procI],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );

            boundaryProcAddressing.set
            (
                procI,
                new labelIOList
                (
                    IOobject
                    (
                        "boundaryRegionAddressing",
                        procMeshes[procI].facesInstance(),
                        procMeshes[procI].meshSubDir,
                        procMeshes[procI],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );
        }
    }

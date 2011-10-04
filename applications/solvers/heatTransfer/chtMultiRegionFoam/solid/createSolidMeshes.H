    forAll(rp.solidRegionNames(), i)
    {
        Info<< "Create solid mesh for region " << rp.solidRegionNames()[i]
            << " for time = " << runTime.timeName() << nl << endl;

        label procI = solidToProc[i];

        procMeshes.set
        (
            procI,
            new fvMesh
            (
                IOobject
                (
                    rp.solidRegionNames()[i],
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

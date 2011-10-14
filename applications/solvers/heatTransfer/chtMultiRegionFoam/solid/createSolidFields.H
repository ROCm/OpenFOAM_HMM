    // Initialise solid field pointer lists
    PtrList<basicSolidThermo> thermoSolid(rp.solidRegionNames().size());
    PtrList<volScalarField> hSolid(rp.solidRegionNames().size());

    // Populate solid field pointer lists
    forAll(rp.solidRegionNames(), i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << rp.solidRegionNames()[i] << nl << endl;

        label procI = solidToProc[i];

        Info<< "    Adding to thermoSolid\n" << endl;
        thermoSolid.set(i, basicSolidThermo::New(procMeshes[procI]));

        if (temperatureCoupled)
        {
            Info<< "    Adding to hSolid\n" << endl;
            hSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "h",
                        runTime.timeName(),
                        procMeshes[procI],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    procMeshes[procI]
                )
            );
        }
    }

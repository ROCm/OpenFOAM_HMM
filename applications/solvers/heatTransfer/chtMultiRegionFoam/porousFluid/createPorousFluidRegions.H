    const wordList porousFluidNames(rp["porousFluid"]);

    PtrList<fvMesh> porousFluidRegions(porousFluidNames.size());

    forAll (porousFluidNames, iPorous)
    {
        const word porousFluidName = porousFluidNames[iPorous];
        Info<< "Create porous fluid region " << porousFluidName
            << nl << endl;

        porousFluidRegions.set
        (
            iPorous,
            new fvMesh
            (
                IOobject
                (
                    porousFluidName,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }

    const wordList porousSolidNames(rp["porousSolid"]);

    PtrList<fvMesh> porousSolidRegions(porousSolidNames.size());

    forAll(porousSolidNames, i)
    {
        Info<< "Create solid mesh for region " << porousSolidNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        porousSolidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    porousSolidNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }

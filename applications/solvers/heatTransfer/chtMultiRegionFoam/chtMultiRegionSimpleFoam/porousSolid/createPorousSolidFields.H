    // Initialise solid field pointer lists
    PtrList<basicSolidThermo> porousSolidThermos(porousSolidRegions.size());
    PtrList<IObasicSourceList> solidHeatSources(porousSolidRegions.size());
    PtrList<volScalarField> betavSolid(porousSolidRegions.size());

    // Populate solid field pointer lists
    forAll(porousSolidRegions, i)
    {
        Info<< "*** Reading porous solid mesh thermophysical "
            << "properties for region "
            << porousSolidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermos\n" << endl;
        porousSolidThermos.set
        (
            i,
            basicSolidThermo::New(porousSolidRegions[i])
        );
        Info<< "    Adding sources\n" << endl;
        solidHeatSources.set
        (
            i,
            new IObasicSourceList(porousSolidRegions[i])
        );

        Info<< "    Adding to betavSolid\n" << endl;
        betavSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "betavSolid",
                    runTime.timeName(),
                    porousSolidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                porousSolidRegions[i]
            )
        );
    }

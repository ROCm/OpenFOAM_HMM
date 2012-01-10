    // Initialise porous field pointer lists
    PtrList<basicRhoThermo> thermoPorous(porousFluidRegions.size());
    PtrList<volScalarField> rhoPorous(porousFluidRegions.size());
    PtrList<volScalarField> kappaPorous(porousFluidRegions.size());
    PtrList<volVectorField> UPorous(porousFluidRegions.size());
    PtrList<surfaceScalarField> phiPorous(porousFluidRegions.size());
    PtrList<compressible::turbulenceModel> turbulencePorous
    (
        porousFluidRegions.size()
    );
    PtrList<volScalarField> pPorous(porousFluidRegions.size());

    List<scalar> initialMassFluidPorous(porousFluidRegions.size());
    List<label> pRefCellFluidPorous(porousFluidRegions.size(),0);
    List<scalar> pRefValueFluidPorous(porousFluidRegions.size(),0.0);

    PtrList<dimensionedScalar> rhoMaxPorous(fluidRegions.size());
    PtrList<dimensionedScalar> rhoMinPorous(fluidRegions.size());

    PtrList<IObasicSourceList> heatPorousSources
    (
        porousFluidRegions.size()
    );

    forAll(porousFluidRegions, i)
    {
        Info<< "Reading fluid mesh thermophysical properties for porous "
            << porousFluidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermoFluid porous\n" << endl;

        thermoPorous.set
        (
            i,
            basicRhoThermo::New(porousFluidRegions[i]).ptr()
        );

        Info<< "    Adding to rhoPorous\n" << endl;
        rhoPorous.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    porousFluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                thermoPorous[i].rho()
            )
        );

        Info<< "    Adding to UPorous\n" << endl;
        UPorous.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    porousFluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                porousFluidRegions[i]
            )
        );

        Info<< "    Adding to phiPorous\n" << endl;
        phiPorous.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    porousFluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(rhoPorous[i]*UPorous[i])
                    & porousFluidRegions[i].Sf()
            )
        );

        Info<< "    Adding turbulence to porous\n" << endl;
        turbulencePorous.set
        (
            i,
            compressible::turbulenceModel::New
            (
                rhoPorous[i],
                UPorous[i],
                phiPorous[i],
                thermoPorous[i]
            ).ptr()
        );

        Info<< "    Adding to kappaFluid\n" << endl;
        kappaPorous.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappaPorous",
                    runTime.timeName(),
                    porousFluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermoPorous[i].Cp()*thermoPorous[i].alpha()
            )
        );

        pPorous.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    porousFluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                porousFluidRegions[i]
            )
        );

        setRefCell
        (
            thermoPorous[i].p(),
            pPorous[i],
            porousFluidRegions[i].solutionDict().subDict("SIMPLE"),
            pRefCellFluidPorous[i],
            pRefValueFluidPorous[i]
        );

        rhoMaxPorous.set
        (
            i,
            new dimensionedScalar
            (
                porousFluidRegions[i].solutionDict().subDict("SIMPLE").lookup
                (
                    "rhoMax"
                )
            )
        );

        rhoMinPorous.set
        (
            i,
            new dimensionedScalar
            (
                porousFluidRegions[i].solutionDict().subDict("SIMPLE").lookup
                (
                    "rhoMin"
                )
            )
        );

        heatPorousSources.set
        (
            i,
            new IObasicSourceList(porousFluidRegions[i])
        );
    }



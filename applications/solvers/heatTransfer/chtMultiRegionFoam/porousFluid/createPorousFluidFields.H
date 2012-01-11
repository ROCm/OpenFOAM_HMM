    // Initialise porous field pointer lists
    PtrList<basicRhoThermo> thermoPorous(porousFluidRegions.size());
    PtrList<volScalarField> rhoPorous(porousFluidRegions.size());
    PtrList<volVectorField> UPorous(porousFluidRegions.size());
    PtrList<surfaceScalarField> phiPorous(porousFluidRegions.size());
    PtrList<volScalarField> KPorous(porousFluidRegions.size());
    PtrList<volScalarField> dpdtPorous(fluidRegions.size());

    PtrList<compressible::turbulenceModel> turbulencePorous
    (
        porousFluidRegions.size()
    );
    PtrList<volScalarField> pPorous(porousFluidRegions.size());

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

        Info<< "    Adding to KPorous\n" << endl;
        KPorous.set
        (
            i,
            new volScalarField
            (
                "KPorous",
                0.5*magSqr(UPorous[i])
            )
        );

        Info<< "    Adding to dpdtPorous\n" << endl;
        dpdtPorous.set
        (
            i,
            new volScalarField
            (
                "dpdtPorous",
                fvc::ddt(thermoPorous[i].p())
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

        heatPorousSources.set
        (
            i,
            new IObasicSourceList(porousFluidRegions[i])
        );
    }



    // Initialise fluid field pointer lists
    PtrList<basicRhoThermo> thermoFluid(rp.fluidRegionNames().size());
    PtrList<volScalarField> rhoFluid(rp.fluidRegionNames().size());
    PtrList<volScalarField> KFluid(rp.fluidRegionNames().size());
    PtrList<volVectorField> UFluid(rp.fluidRegionNames().size());
    PtrList<surfaceScalarField> phiFluid(rp.fluidRegionNames().size());
    PtrList<uniformDimensionedVectorField> gFluid(rp.fluidRegionNames().size());
    PtrList<compressible::turbulenceModel> turbulence
    (
        rp.fluidRegionNames().size()
    );
    PtrList<volScalarField> p_rghFluid(rp.fluidRegionNames().size());
    PtrList<volScalarField> ghFluid(rp.fluidRegionNames().size());
    PtrList<surfaceScalarField> ghfFluid(rp.fluidRegionNames().size());
    PtrList<radiation::radiationModel> radiation(rp.fluidRegionNames().size());
    PtrList<volScalarField> DpDtFluid(rp.fluidRegionNames().size());

    List<scalar> initialMassFluid(rp.fluidRegionNames().size());

    // Populate fluid field pointer lists
    forAll(rp.fluidRegionNames(), i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << rp.fluidRegionNames()[i] << nl << endl;

        label procI = fluidToProc[i];


        Info<< "    Adding to thermoFluid\n" << endl;
        thermoFluid.set
        (
            i,
            basicRhoThermo::New(procMeshes[procI]).ptr()
        );

        Info<< "    Adding to rhoFluid\n" << endl;
        rhoFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    procMeshes[procI],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                thermoFluid[i].rho()
            )
        );

        Info<< "    Adding to KFluid\n" << endl;
        KFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "K",
                    runTime.timeName(),
                    procMeshes[procI],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermoFluid[i].Cp()*thermoFluid[i].alpha()
            )
        );

        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    procMeshes[procI],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                procMeshes[procI]
            )
        );

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    procMeshes[procI],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(rhoFluid[i]*UFluid[i])
                    & procMeshes[procI].Sf()
            )
        );

        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    procMeshes[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        Info<< "    Adding to turbulence\n" << endl;
        turbulence.set
        (
            i,
            autoPtr<compressible::turbulenceModel>
            (
                compressible::turbulenceModel::New
                (
                    rhoFluid[i],
                    UFluid[i],
                    phiFluid[i],
                    thermoFluid[i]
                )
            ).ptr()
        );

        Info<< "    Adding to ghFluid\n" << endl;
        ghFluid.set
        (
            i,
            new volScalarField("gh", gFluid[i] & procMeshes[procI].C())
        );

        Info<< "    Adding to ghfFluid\n" << endl;
        ghfFluid.set
        (
            i,
            new surfaceScalarField("ghf", gFluid[i] & procMeshes[procI].Cf())
        );

        p_rghFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p_rgh",
                    runTime.timeName(),
                    procMeshes[procI],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                procMeshes[procI]
            )
        );

        // Force p_rgh to be consistent with p
        p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

        radiation.set
        (
            i,
            radiation::radiationModel::New(thermoFluid[i].T())
        );

        initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

        Info<< "    Adding to DpDtFluid\n" << endl;
        DpDtFluid.set
        (
            i,
            new volScalarField
            (
                "DpDt",
                fvc::DDt
                (
                    surfaceScalarField
                    (
                        "phiU",
                        phiFluid[i]/fvc::interpolate(rhoFluid[i])
                    ),
                    thermoFluid[i].p()
                )
            )
        );
    }

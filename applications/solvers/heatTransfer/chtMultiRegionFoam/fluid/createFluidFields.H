    // Initialise fluid field pointer lists
    PtrList<basicThermo> thermoFluid(fluidRegions.size());
    PtrList<volScalarField> rhoFluid(fluidRegions.size());
    PtrList<volScalarField> KFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
    PtrList<compressible::turbulenceModel> turbulence(fluidRegions.size());
    PtrList<volScalarField> DpDtFluid(fluidRegions.size());
    PtrList<volScalarField> ghFluid(fluidRegions.size());
    PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
    PtrList<volScalarField> pdFluid(fluidRegions.size());

    List<scalar> initialMassFluid(fluidRegions.size());

    dimensionedScalar pRef
    (
        "pRef",
        dimensionSet(1, -1, -2, 0, 0),
        rp.lookup("pRef")
    );

    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to pdFluid\n" << endl;
        pdFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "pd",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to thermoFluid\n" << endl;
        thermoFluid.set
        (
            i,
            basicThermo::New(fluidRegions[i]).ptr()
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
                    fluidRegions[i],
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
                    fluidRegions[i],
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
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
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
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(rhoFluid[i]*UFluid[i])
                    & fluidRegions[i].Sf()
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

        Info<< "    Adding to DpDtFluid\n" << endl;
        DpDtFluid.set
        (
            i,
            new volScalarField
            (
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

        const dictionary& environmentalProperties =
            fluidRegions[i].lookupObject<IOdictionary>
                ("environmentalProperties");
        dimensionedVector g(environmentalProperties.lookup("g"));

        Info<< "    Adding to ghFluid\n" << endl;
        ghFluid.set
        (
            i,
            new volScalarField
            (
                "gh",
                 g & fluidRegions[i].C()
            )
        );
        ghfFluid.set
        (
            i,
            new surfaceScalarField
            (
                "ghf",
                 g & fluidRegions[i].Cf()
            )
        );

        Info<< "    Updating p from pd\n" << endl;
        thermoFluid[i].p() == pdFluid[i] + rhoFluid[i]*ghFluid[i] + pRef;
        thermoFluid[i].correct();

        initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();
    }



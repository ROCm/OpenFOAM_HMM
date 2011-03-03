    Info<< "\nReading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rhoInfValue
    (
        transportProperties.lookup("rhoInf")
    );

    volScalarField rhoInf
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        rhoInfValue
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    Info<< "Creating turbulence model\n" << endl;

    singlePhaseTransportModel laminarTransport(U, phi);

    const volScalarField nu(laminarTransport.nu());

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        nu*rhoInfValue
    );

    word kinematicCloudName("kinematicCloud");
    args.optionReadIfPresent("cloudName", kinematicCloudName);

    Info<< "Constructing kinematicCloud " << kinematicCloudName << endl;
    basicKinematicCollidingCloud kinematicCloud
    (
        kinematicCloudName,
        rhoInf,
        U,
        mu,
        g
    );

    IOobject Hheader
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

    autoPtr<volVectorField> HPtr_;

    if (Hheader.headerOk())
    {
        Info<< "\nReading field H\n" << endl;

        HPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "H",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

    IOobject HdotGradHheader
    (
        "HdotGradH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

    autoPtr<volVectorField> HdotGradHPtr_;

    if (HdotGradHheader.headerOk())
    {
        Info<< "Reading field HdotGradH" << endl;

        HdotGradHPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "HdotGradH",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

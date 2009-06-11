    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicThermo> thermo
    (
        basicThermo::New(mesh)
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo->rho()
    );

    Info<< "\nReading field U\n" << endl;
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

    #include "compressibleCreatePhi.H"

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo()
        )
    );

    word kinematicCloudName("kinematicCloud");
    args.optionReadIfPresent("cloudName", kinematicCloudName);

    Info<< "Constructing kinematicCloud " << kinematicCloudName << endl;
    basicKinematicCloud kinematicCloud
    (
        kinematicCloudName,
        rho,
        U,
        thermo().mu(),
        g
    );

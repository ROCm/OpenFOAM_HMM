    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
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

    // Create p before the thermo
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh
    );

    // Creating e based thermo
    autoPtr<twoPhaseMixtureEThermo> thermo;
    thermo.set(new twoPhaseMixtureEThermo(U, phi));

    // Create mixture and
    Info<< "Creating temperaturePhaseChangeTwoPhaseMixture\n" << endl;
    autoPtr<temperaturePhaseChangeTwoPhaseMixture> mixture =
        temperaturePhaseChangeTwoPhaseMixture::New(thermo(), mesh);


    volScalarField& T = thermo->T();
    volScalarField& e = thermo->he();

    // Correct e from T and alpha
    thermo->correct();

    volScalarField& alpha1(thermo->alpha1());
    volScalarField& alpha2(thermo->alpha2());

    const dimensionedScalar& rho1 = thermo->rho1();
    const dimensionedScalar& rho2 = thermo->rho2();

    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + alpha2*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


    // Construct interface from alpha1 distribution
    interfaceProperties interface
    (
        alpha1,
        U,
        thermo->transportPropertiesDict()
    );

    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, thermo())
    );


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    //Update p with rho
    p = p_rgh + rho*gh;

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    // Turbulent Prandtl number
    dimensionedScalar Prt("Prt", dimless, thermo->transportPropertiesDict());

    volScalarField kappaEff
    (
        IOobject
        (
            "kappaEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo->kappa()
    );

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));

    Info<< "Creating field pDivU\n" << endl;
    volScalarField pDivU
    (
        IOobject
        (
            "pDivU",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("pDivU", p.dimensions()/dimTime, 0)
    );

    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.rho()
    );

    volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();
    const volScalarField& psi = thermo.psi();

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

    #include "compressibleCreatePhi.H"


    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    p = p_rgh + rho*gh;
    thermo.correct();
    rho = thermo.rho();
    p_rgh = p - rho*gh;

    label p_rghRefCell = 0;
    scalar p_rghRefValue = 0.0;
    setRefCell
    (
        p_rgh,
        mesh.solutionDict().subDict("SIMPLE"),
        p_rghRefCell,
        p_rghRefValue
    );

    scalar pRefValue = 0.0;

    if (p_rgh.needReference())
    {
        pRefValue = readScalar
        (
            mesh.solutionDict().subDict("SIMPLE").lookup("pRefValue")
        );

        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, p_rghRefCell)
        );
    }

    dimensionedScalar initialMass = fvc::domainIntegrate(rho);

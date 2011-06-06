    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<psiChemistryModel> pChemistry
    (
        psiChemistryModel::New(mesh)
    );
    psiChemistryModel& chemistry = pChemistry();

    hsCombustionThermo& thermo = chemistry.thermo();

    SLGThermo slgThermo(mesh, thermo);

    basicMultiComponentMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    const word inertSpecie(thermo.lookup("inertSpecie"));

    Info<< "Creating field rho\n" << endl;
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
        thermo.rho()
    );

    volScalarField& p = thermo.p();
    volScalarField& hs = thermo.hs();
    const volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();

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

    Info<< "Creating field kappa\n" << endl;
    DimensionedField<scalar, volMesh> kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Creating field DpDt\n" << endl;
    volScalarField DpDt
    (
        "DpDt",
        fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p)
    );


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());

    surfaceScalarField ghf("ghf", g & mesh.Cf());

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

    // Force p_rgh to be consistent with p
    p_rgh = p - rho*gh;

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(hs);

    IOdictionary additionalControlsDict
    (
        IOobject
        (
            "additionalControls",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Switch solvePrimaryRegion
    (
        additionalControlsDict.lookup("solvePrimaryRegion")
    );

    DimensionedField<scalar, volMesh> chemistrySh
    (
        IOobject
        (
            "chemistry::Sh",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0)
    );

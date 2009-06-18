    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<hCombustionThermo> thermo
    (
        hCombustionThermo::New(mesh)
    );

    combustionMixture& composition = thermo->composition();
    PtrList<volScalarField>& Y = composition.Y();

    word inertSpecie(thermo->lookup("inertSpecie"));

    volScalarField& p = thermo->p();
    volScalarField& h = thermo->h();
    const volScalarField& T = thermo->T();
    const volScalarField& psi = thermo->psi();

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

    // lagrangian effective density field - used externally (optional)
    volScalarField rhoEffLagrangian
    (
        IOobject
        (
            "rhoEffLagrangian",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimDensity, 0.0)
    );

    // dynamic pressure field - used externally (optional)
    volScalarField pDyn
    (
        IOobject
        (
            "pDyn",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
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

#   include "compressibleCreatePhi.H"

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
            thermo()
        )
    );

    Info<< "Creating field DpDt\n" << endl;
    volScalarField DpDt =
        fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    Info<< "\nConstructing explicit enthalpy cell source" << endl;
    timeActivatedExplicitCellSource enthalpySource
    (
        "enthalpySource",
        mesh,
        h.dimensions()*phi.dimensions()/mesh.V().dimensions()
    );

    Info << "Constructing chemical mechanism" << endl;
    chemistryModel chemistry
    (
        thermo(),
        rho
    );

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll (Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(h);

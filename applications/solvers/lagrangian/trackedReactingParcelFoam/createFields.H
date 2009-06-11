    Info<< "Reading thermophysical properties" << nl << endl;

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

    Info<< "Reading field U" << nl << endl;
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

    Info<< "Creating turbulence model" << nl << endl;
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

    Info<< "Creating field DpDt" << nl << endl;
    volScalarField DpDt =
        fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    Info << "Constructing chemical mechanism" << nl << endl;
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

    Info<< "Creating porous zones" << nl << endl;
    porousZones pZones(mesh);
    Switch pressureImplicitPorosity(false);

    label nUCorr = 0;
    if (pZones.size())
    {
        // nUCorrectors for pressureImplicitPorosity
        if (mesh.solutionDict().subDict("PISO").found("nUCorrectors"))
        {
            mesh.solutionDict().subDict("PISO").lookup("nUCorrectors")
                >> nUCorr;
        }

        if (nUCorr > 0)
        {
            pressureImplicitPorosity = true;
        }
    }


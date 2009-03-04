
    Info<< nl << "Reading field U" << endl;
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

    Info<< nl << "Reading field T" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field rhoN" << endl;
    volScalarField rhoN
    (
        IOobject
        (
            "rhoN",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field rhoM" << endl;
    volScalarField rhoM
    (
        IOobject
        (
            "rhoM",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field rhoNdsmc (dsmc particle density)" << endl;
    volScalarField dsmcRhoN
    (
        IOobject
        (
            "dsmcRhoN",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl
        << "Fields required to record solid surface forces and heat flux:"
        << endl;

    Info<< nl << "Reading field q (surface heat transfer)" << endl;
    volScalarField q
    (
        IOobject
        (
            "q",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field fD (surface force density" << endl;
    volVectorField fD
    (
        IOobject
        (
            "fD",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Constructing dsmcCloud " << endl;

dsmcCloud dsmc("dsmc", T, U);
    //, rhoN, rhoM, dsmcRhoN, q, f);

scalar timeToStartAveraging = runTime.value();

volVectorField Umean
(
    IOobject
    (
        "Umean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U
);

volScalarField rhoMean
(
    IOobject
    (
        "rhoMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho
);

volScalarField pMean
(
    IOobject
    (
        "pMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p
);

volScalarField gammaMean
(
    IOobject
    (
        "gammaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gamma
);

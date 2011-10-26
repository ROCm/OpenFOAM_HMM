Info<< "Reading thermal properties\n" << endl;

IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Switch thermalStress(thermalProperties.lookup("thermalStress"));

volScalarField threeKalpha
(
    IOobject
    (
        "threeKalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 2, -2 , -1, 0), 0.0)
);


volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 2, -1 , 0, 0), 0.0)
);


if (thermalStress)
{
    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(0),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 2, -2 , -1, 0), 0.0)
    );

    const dictionary& CDict(thermalProperties.subDict("C"));
    word CType(CDict.lookup("C"));
    if (CType == "CInf")
    {
        C = CDict.lookup("CInf");
    }


    volScalarField rhoK
    (
        IOobject
        (
            "k",
            runTime.timeName(0),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1, 1, -3 , -1, 0), 0.0)
    );

    const dictionary& kDict(thermalProperties.subDict("k"));
    word kType(kDict.lookup("k"));
    if (kType == "kInf")
    {
        rhoK = kDict.lookup("kInf");
    }

    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(0),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 0, 0 , -1, 0), 0.0)
    );

    const dictionary& alphaDict(thermalProperties.subDict("alpha"));
    word alphaType(alphaDict.lookup("alpha"));

    if (alphaType == "alphaInf")
    {
        alpha = alphaDict.lookup("alphaInf");
    }

    Info<< "Normalising k : k/rho\n" << endl;
    volScalarField k = rhoK/rho;

    Info<< "Calculating thermal coefficients\n" << endl;

    threeKalpha = threeK*alpha;
    DT = k/C;

}

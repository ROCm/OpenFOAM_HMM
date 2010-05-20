    Info<< "Reading thermodynamicProperties\n" << endl;

    IOdictionary thermodynamicProperties
    (
        IOobject
        (
            "thermodynamicProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar R
    (
        thermodynamicProperties.lookup("R")
    );

    dimensionedScalar Cv
    (
        thermodynamicProperties.lookup("Cv")
    );

    dimensionedScalar Cp = Cv + R;

    dimensionedScalar gamma = Cp/Cv;

    dimensionedScalar Pr
    (
        dimensionedScalar::lookupOrDefault
        (
            "Pr",
            thermodynamicProperties,
            1.0
        )
    );


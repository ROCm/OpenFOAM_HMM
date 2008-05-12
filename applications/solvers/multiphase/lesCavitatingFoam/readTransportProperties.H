    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar mul
    (
        transportProperties.lookup("mul")
    );

    dimensionedScalar muv
    (
        transportProperties.lookup("muv")
    );

    IOdictionary ppProperties
    (
        IOobject
        (
            "ppProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar preAlphaExp
    (
        readScalar(ppProperties.lookup("preAlphaExp"))
    );

    scalar alpha1Max
    (
        readScalar(ppProperties.lookup("alpha1Max"))
    );

    scalar expMax
    (
        readScalar(ppProperties.lookup("expMax"))
    );

    dimensionedScalar g0
    (
        ppProperties.lookup("g0")
    );

    Switch packingLimiter
    (
        ppProperties.lookup("packingLimiter")
    );

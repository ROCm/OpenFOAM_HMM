{
    volScalarField rhoEff
    (
        IOobject
        (
            "rhoEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1.0 - beta*(T - TRef)
    );
    rhoEff.write();

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pd + rhoEff*(g & mesh.C()) + pRef
    );
    p.write();
}

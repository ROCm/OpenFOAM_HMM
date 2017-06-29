#if defined(MOVING_MESH)
    Info<< "Creating field dpdt for moving meshes\n" << endl;

    // Note: set to read if present and auto write to simplify dpdt correction
    // by meshPhi
    volScalarField dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt(p)
    );
#else
    Info<< "Creating field dpdt\n" << endl;

    volScalarField dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::ddt(p)
    );
#endif

if (!thermo.dpdt())
{
    dpdt == dimensionedScalar("0", dpdt.dimensions(), 0);
    dpdt.writeOpt() = IOobject::NO_WRITE;
}

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

if (!thermo.dpdt())
{
    dpdt == dimensionedScalar("0", dpdt.dimensions(), 0);
}

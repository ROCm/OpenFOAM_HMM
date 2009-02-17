Info<< "\nReading environmentalProperties" << endl;

IOdictionary environmentalProperties
(
    IOobject
    (
        "environmentalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedVector g(environmentalProperties.lookup("g"));
const Switch rotating(environmentalProperties.lookup("rotating"));
const dimensionedVector Omega =
    rotating ? environmentalProperties.lookup("Omega")
             : dimensionedVector("Omega", -dimTime, vector(0,0,0));
const dimensionedScalar magg = mag(g);
const dimensionedVector gHat = g/magg;

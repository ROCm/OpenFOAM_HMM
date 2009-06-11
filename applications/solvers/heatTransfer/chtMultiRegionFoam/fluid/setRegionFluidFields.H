    const fvMesh& mesh = fluidRegions[i];

    basicThermo& thermo = thermoFluid[i];
    volScalarField& rho = rhoFluid[i];
    volScalarField& K = KFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    compressible::turbulenceModel& turb = turbulence[i];
    volScalarField& DpDt = DpDtFluid[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();
    volScalarField& h = thermo.h();

    const dimensionedScalar massIni("massIni", dimMass, initialMassFluid[i]);

    const dictionary& environmentalProperties =
        fluidRegions[i].lookupObject<IOdictionary>
        ("environmentalProperties");

    const dimensionedVector g(environmentalProperties.lookup("g"));

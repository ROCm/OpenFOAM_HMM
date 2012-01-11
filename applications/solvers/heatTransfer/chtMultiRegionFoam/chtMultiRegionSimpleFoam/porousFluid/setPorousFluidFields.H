
    const fvMesh& porousMesh = porousFluidRegions[i];

    basicRhoThermo& porousThermo = thermoPorous[i];
    volScalarField& porousRho = rhoPorous[i];

    volVectorField& porousU = UPorous[i];
    surfaceScalarField& porousPhi = phiPorous[i];

    compressible::turbulenceModel& turbPorous = turbulencePorous[i];

    volScalarField& porousP = porousThermo.p();
    const volScalarField& porousPsi = porousThermo.psi();
    volScalarField& porousH = porousThermo.h();

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluidPorous[i]
    );

    IObasicSourceList& porousSources = heatPorousSources[i];

    const label pRefCell = pRefCellFluidPorous[i];
    const scalar pRefValue = pRefValueFluidPorous[i];

    const scalar rhoMax = rhoMaxPorous[i].value();
    const scalar rhoMin = rhoMinPorous[i].value();



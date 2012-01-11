
    fvMesh& porousMesh = porousFluidRegions[i];

    basicRhoThermo& porousThermo = thermoPorous[i];
    volScalarField& porousRho = rhoPorous[i];

    volVectorField& porousU = UPorous[i];
    surfaceScalarField& porousPhi = phiPorous[i];

    compressible::turbulenceModel& turbPorous = turbulencePorous[i];
    volScalarField& porousK = KPorous[i];
    volScalarField& porousdpdt = dpdtPorous[i];

    volScalarField& porousP = porousThermo.p();
    const volScalarField& porousPsi = porousThermo.psi();
    volScalarField& porousH = porousThermo.h();

    IObasicSourceList& porousSources = heatPorousSources[i];

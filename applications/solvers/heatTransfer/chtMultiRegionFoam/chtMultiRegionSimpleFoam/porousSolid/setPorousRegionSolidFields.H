    const fvMesh& mesh = porousSolidRegions[i];

    basicSolidThermo& thermo = porousSolidThermos[i];
    const volScalarField& betav = betavSolid[i];

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> tkappa = betav*thermo.K();
    //tmp<volSymmTensorField> tkappa = thermo.directionalK()*betav;

    const volScalarField& kappa = tkappa();
    //const volSymmTensorField& K = tK();

    volScalarField& T = thermo.T();

    IObasicSourceList& sources = solidHeatSources[i];

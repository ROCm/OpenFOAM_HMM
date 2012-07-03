    fvMesh& mesh = porousSolidRegions[i];

    solidThermo& thermo = porousSolidThermos[i];
    const volScalarField& betav = betavSolid[i];

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> tkappa = thermo.kappa();
    //tmp<volSymmTensorField> tkappa = thermo.directionalKappa()*betav;

    const volScalarField& kappa = tkappa();
    //const volSymmTensorField& K = tK();

    //tmp<volScalarField> trhoCp = cp*rho;
    //const volScalarField& rhoCp = trhoCp();

    tmp<volScalarField> talpha = thermo.alpha();
    const volScalarField& alpha = talpha();

    volScalarField& h = thermo.he();

    IObasicSourceList& sources = solidHeatSources[i];

    fvMesh& mesh = solidRegions[i];
    solidThermo& thermo = thermos[i];

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> talpha = thermo.alpha();
    const volScalarField& alpha = talpha();
    tmp<volScalarField> tkappa = thermo.kappa();
    const volScalarField& kappa = tkappa();

    //tmp<volSymmTensorField> tkappa = thermo.directionalKappa();
    //const volSymmTensorField& kappa = tkappa();

    volScalarField& h = thermo.he();

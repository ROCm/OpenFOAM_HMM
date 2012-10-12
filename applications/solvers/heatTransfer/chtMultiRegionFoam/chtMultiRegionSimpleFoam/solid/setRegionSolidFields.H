    fvMesh& mesh = solidRegions[i];
    solidThermo& thermo = thermos[i];
    const radiation::radiationModel& radiation = radiations[i];

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> tkappa = thermo.kappa();
    //tmp<volSymmTensorField> tkappa = thermo.directionalkappa();
    const volScalarField& kappa = tkappa();

    tmp<volScalarField> talpha = thermo.alpha();
    const volScalarField& alpha = talpha();

    volScalarField& h = thermo.he();

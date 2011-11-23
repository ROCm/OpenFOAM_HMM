    fvMesh& mesh = solidRegions[i];
    basicSolidThermo& thermo = thermos[i];

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> tK = thermo.K();
    //tmp<volSymmTensorField> tK = thermo.directionalK();
    const volScalarField& K = tK();
    //const volSymmTensorField& K = tK();

    volScalarField& T = thermo.T();

{
    volScalarField kappaEff
    (
        "kappaEff",
        turbulence->nu() + turbulence->nut()/Pr
    );

    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::Sp(fvc::div(phi), T)
      - fvm::laplacian(kappaEff, T)
    );

    TEqn.relax();

    eqnResidual = TEqn.solve().initialResidual();
    maxResidual = max(eqnResidual, maxResidual);
}

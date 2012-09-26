{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ep", p/rho))
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
    );

    EEqn.relax();
    EEqn.solve();

    thermo.correct();
}

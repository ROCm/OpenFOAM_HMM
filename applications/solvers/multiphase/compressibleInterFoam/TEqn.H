{
    solve
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvm::laplacian(thermo.alphaEff(rho*turbulence->nut()), T)
      + (
            p*fvc::div(phi)
          + fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
        )*(alpha1/thermo.thermo1().Cv() + alpha2/thermo.thermo2().Cv())
    );

    thermo.correct();
}

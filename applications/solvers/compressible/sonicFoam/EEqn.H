{
    solve
    (
        fvm::ddt(rho, e) + fvm::div(phi, e)
      + fvc::ddt(rho, K) + fvc::div(phi, volScalarField("Ekp", K + p/rho))
      - fvm::laplacian(turbulence->alphaEff(), e)
    );

    thermo.correct();
}

{
    solve
    (
        fvm::ddt(rho, e)
      + fvm::div(phi, e)
      - fvm::laplacian(turbulence->alphaEff(), e)
     ==
      - (fvc::ddt(rho, K) + fvc::div(phi, volScalarField("Ekp", K + p/rho)))
    );

    thermo.correct();
}

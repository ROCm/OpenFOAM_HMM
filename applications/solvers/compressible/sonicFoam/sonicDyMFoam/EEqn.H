{
    solve
    (
        fvm::ddt(rho, e) + fvm::div(phi, e)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + fvc::div(phi/fvc::interpolate(rho) + mesh.phi(), p, "div(phiv,p)")
      - fvm::laplacian(turbulence->alphaEff(), e)
    );

    thermo.correct();
}

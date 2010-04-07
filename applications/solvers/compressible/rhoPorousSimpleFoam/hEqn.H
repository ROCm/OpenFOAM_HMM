{
    fvScalarMatrix hEqn
    (
        fvm::div(phi, h)
      - fvm::Sp(fvc::div(phi), h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        fvc::div(phi/fvc::interpolate(rho)*fvc::interpolate(p, "div(U,p)"))
      - p*fvc::div(phi/fvc::interpolate(rho))
    );

    pZones.addEnthalpySource(thermo, rho, hEqn);

    hEqn.relax();

    hEqn.solve();

    thermo.correct();
}

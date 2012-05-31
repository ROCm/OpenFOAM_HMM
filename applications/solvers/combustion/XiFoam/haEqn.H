{
    fvScalarMatrix haEqn
    (
        fvm::ddt(rho, ha)
      + mvConvection->fvmDiv(phi, ha)
      - fvm::laplacian(turbulence->alphaEff(), ha)
     ==
        dpdt
      - (fvc::ddt(rho, K) + fvc::div(phi, K))
    );

    haEqn.relax();
    haEqn.solve();

    thermo.correct();
}

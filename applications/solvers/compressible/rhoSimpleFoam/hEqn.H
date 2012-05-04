{
    volScalarField K("K", 0.5*magSqr(U));

    fvScalarMatrix hEqn
    (
        fvm::div(phi, h)
      - fvm::Sp(fvc::div(phi), h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        fvc::div(phi)*K - fvc::div(phi, K)
    );

    hEqn.relax();
    hEqn.solve();

    thermo.correct();
}

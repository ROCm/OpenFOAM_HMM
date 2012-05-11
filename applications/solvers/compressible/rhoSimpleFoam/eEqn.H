{
    // Kinetic + pressure energy
    volScalarField Ekp("Ekp", 0.5*magSqr(U) + p/rho);

    fvScalarMatrix eEqn
    (
        fvm::div(phi, e)
      - fvm::Sp(fvc::div(phi), e)
      - fvm::laplacian(turbulence->alphaEff(), e)
     ==
        fvc::div(phi)*Ekp - fvc::div(phi, Ekp)
    );

    eEqn.relax();
    eEqn.solve();

    thermo.correct();
}

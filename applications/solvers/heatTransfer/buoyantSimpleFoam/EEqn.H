{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        radiation->Sh(thermo)
      + sources(rho, he)
    );

    EEqn.relax();

    sources.constrain(EEqn);

    EEqn.solve();

    thermo.correct();
    radiation->correct();
}

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
      - fvm::laplacian(turb.alphaEff(), he)
     ==
        rad.Sh(thermo)
      + sources(rho, he)
    );

    EEqn.relax();

    sources.constrain(EEqn);

    EEqn.solve();

    sources.correct(he);

    thermo.correct();
    rad.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}

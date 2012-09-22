{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ep", p/rho))
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        parcels.Sh(he)
      + surfaceFilm.Sh()
      + radiation->Sh(thermo)
      + combustion->Sh()
    );

    EEqn.relax();
    EEqn.solve();

    thermo.correct();
    radiation->correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}

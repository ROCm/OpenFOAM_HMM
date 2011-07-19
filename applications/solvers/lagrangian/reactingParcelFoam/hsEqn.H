{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
        DpDt
     +  parcels.Sh(hs)
     +  radiation->Shs(thermo)
     +  combustion->Sh()
    );

    hEqn.relax();

    hEqn.solve();

    thermo.correct();

    radiation->correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}

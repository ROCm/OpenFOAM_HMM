    // Relative momentum predictor

    tmp<fvVectorMatrix> UrelEqn
    (
        fvm::div(phi, Urel)
      + turbulence->divDevReff(Urel)
      + SRF->Su()
    );

    sources.apply(UrelEqn());

    UrelEqn().relax();

    solve(UrelEqn() == -fvc::grad(p));

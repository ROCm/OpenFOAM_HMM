    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    sources.apply(UEqn());

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));

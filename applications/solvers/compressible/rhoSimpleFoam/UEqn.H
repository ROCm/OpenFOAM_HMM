    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        sources(rho, U)
    );

    UEqn().relax();

    sources.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    sources.correct(U);

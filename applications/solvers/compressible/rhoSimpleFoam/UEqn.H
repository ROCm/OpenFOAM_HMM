    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        sources(rho, U)
    );

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));

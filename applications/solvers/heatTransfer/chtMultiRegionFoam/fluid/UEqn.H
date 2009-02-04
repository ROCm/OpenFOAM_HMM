    // Solve the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turb.divDevRhoReff(U)
    );

    UEqn().relax();

    if (momentumPredictor)
    {
        solve(UEqn() == -fvc::grad(pd) - fvc::grad(rho)*gh);
    }

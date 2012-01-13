    // Solve the Momentum equation
    tmp<fvVectorMatrix> porousUEqn
    (
        fvm::div(porousPhi, porousU)
      + turbPorous.divDevRhoReff(porousU)
      + porousSources(porousRho, porousU)
    );

    porousUEqn().relax();

    solve(porousUEqn() == -fvc::grad(porousP));

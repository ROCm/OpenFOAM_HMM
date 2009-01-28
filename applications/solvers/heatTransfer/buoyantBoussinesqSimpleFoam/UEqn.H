    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      - fvm::Sp(fvc::div(phi), U)
      + turbulence->divDevReff(U)
    );

    UEqn().relax();

    eqnResidual = solve
    (
        UEqn()
      ==
       -fvc::grad(pd)
      + beta*gh*fvc::grad(T)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);

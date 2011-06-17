    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    // Add resistance on the actuation disks
    actuationDisks.addSu(UEqn());

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));

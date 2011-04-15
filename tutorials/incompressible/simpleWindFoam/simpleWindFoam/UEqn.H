    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    UEqn().relax();

    // Add resistance on the actuation disks
    actuationDisks.addSu(UEqn());

    solve(UEqn() == -fvc::grad(p));

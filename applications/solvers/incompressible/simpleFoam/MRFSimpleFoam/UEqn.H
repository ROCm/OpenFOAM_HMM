    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
        + turbulence->divDevReff(U)
    );
    mrfZones.addCoriolis(UEqn());

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));

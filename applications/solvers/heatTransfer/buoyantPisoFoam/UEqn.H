    // Solve the Momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
    );

    UEqn.relax();

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         ==
           -fvc::reconstruct
            (
                (
                    fvc::snGrad(pd)
                  + ghf*fvc::snGrad(rho)
                ) * mesh.magSf()
            )
        );
    }

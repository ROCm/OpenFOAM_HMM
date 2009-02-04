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
         solve
         (
            UEqn()
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

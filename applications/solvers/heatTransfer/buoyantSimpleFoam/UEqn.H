    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      - fvm::Sp(fvc::div(phi), U)
      + turbulence->divDevRhoReff(U)
    );

    UEqn().relax();

    eqnResidual = solve
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
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);


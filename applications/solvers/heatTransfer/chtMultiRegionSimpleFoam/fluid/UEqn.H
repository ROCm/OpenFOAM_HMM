    // Solve the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turb.divDevRhoReff(U)
    );

    UEqn().relax();

    solve
    (
        UEqn()
     ==
        fvc::reconstruct
        (
            fvc::interpolate(rho)*(g & mesh.Sf())
          - fvc::snGrad(p)*mesh.magSf()
        )
    );

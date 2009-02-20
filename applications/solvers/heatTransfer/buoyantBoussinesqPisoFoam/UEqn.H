    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    UEqn().relax();

    solve
    (
        UEqn()
      ==
       -fvc::reconstruct
        (
            (
                fvc::snGrad(pd)
              - betaghf*fvc::snGrad(T)
            ) * mesh.magSf()
        )
    );

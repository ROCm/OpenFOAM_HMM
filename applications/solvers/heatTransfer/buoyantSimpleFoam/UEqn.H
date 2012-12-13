    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        sources(rho, U)
    );

    UEqn().relax();

    sources.constrain(UEqn());

    if (simple.momentumPredictor())
    {
        solve
        (
            UEqn()
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        sources.correct(U);
    }

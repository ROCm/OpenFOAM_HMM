    surfaceScalarField muf =
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut());

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muf, U)
      - (fvc::grad(U) & fvc::grad(muf))
    //- fvc::div(muf*(mesh.Sf() & fvc::interpolate(fvc::grad(U)().T())))
    );

    UEqn.relax();

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(pd)
                ) * mesh.magSf()
            )
        );
    }

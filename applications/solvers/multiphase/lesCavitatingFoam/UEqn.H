    surfaceScalarField gammaf = fvc::interpolate(gamma);
    surfaceScalarField muf
        (
            "muf",
            gammaf*muv + (1.0 - gammaf)*mul
          + fvc::interpolate(rho*turbulence->nuSgs())
        );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      - fvm::laplacian(muf, U)
    //- (fvc::grad(U) & fvc::grad(muf))
      - fvc::div(muf*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        thermoCloud1.SU(U)
      + kinematicCloud1.SU(U)
      + rho.dimensionedInternalField()*g
    );

    UEqn.relax();

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }

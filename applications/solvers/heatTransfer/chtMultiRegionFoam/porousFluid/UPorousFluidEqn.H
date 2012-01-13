
// Solve the Momentum equation
    tmp<fvVectorMatrix> porousUEqn
    (
        fvm::ddt(porousRho, porousU)
      + fvm::div(porousPhi, porousU)
      + turbPorous.divDevRhoReff(porousU)
      + porousSources(porousRho, porousU)
    );

    porousUEqn().relax();

    volScalarField rAUPorous(1.0/porousUEqn().A());

    if (momentumPredictor)
    {
        solve(porousUEqn() == -fvc::grad(porousP));
    }
    else
    {
        porousU = rAUPorous*(porousUEqn().H() - fvc::grad(porousP));
        porousU.correctBoundaryConditions();
    }


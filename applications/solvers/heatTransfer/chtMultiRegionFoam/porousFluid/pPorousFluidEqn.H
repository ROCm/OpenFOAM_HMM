porousRho = porousThermo.rho();

porousU = rAUPorous*porousUEqn().H();

if (nCorr <= 1)
{
    porousUEqn.clear();
}

porousPhi =
    fvc::interpolate(porousRho)*
    (
        (fvc::interpolate(porousU) & porousMesh.Sf())
        + fvc::ddtPhiCorr(rAUPorous, porousRho, porousU, porousPhi)
    );

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::ddt(porousPsi, porousP)
      + fvc::div(porousPhi)
      - fvm::laplacian(porousRho*rAUPorous, porousP)
    );

    pEqn.solve
    (
        porousMesh.solver
        (
            porousP.select
            (
                oCorr == nOuterCorr-1
                && corr == nCorr-1
                && nonOrth == nNonOrthCorr
            )
        )
    );

    if (nonOrth == nNonOrthCorr)
    {
        porousPhi += pEqn.flux();
    }
}

solve(fvm::ddt(porousRho) + fvc::div(porousPhi));

// Explicitly relax pressure for momentum corrector
porousP.relax();

// Recalculate density from the relaxed pressure
porousRho = porousThermo.rho();


porousU -= rAUPorous*fvc::grad(porousP);
porousU.correctBoundaryConditions();
porousK = 0.5*magSqr(porousU);

// Update pressure time derivative
porousdpdt = fvc::ddt(porousP);


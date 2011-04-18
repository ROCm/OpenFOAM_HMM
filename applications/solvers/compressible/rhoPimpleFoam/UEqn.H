// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
);

UEqn().relax();

volScalarField rAU(1.0/UEqn().A());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));
}
else
{
    U = rAU*(UEqn().H() - fvc::grad(p));
    U.correctBoundaryConditions();
}

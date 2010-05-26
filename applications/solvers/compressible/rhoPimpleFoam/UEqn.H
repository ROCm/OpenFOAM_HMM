// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
);

UEqn().relax();

volScalarField rUA = 1.0/UEqn().A();

if (momentumPredictor)
{
    solve(UEqn() == -fvc::grad(p));
}
else
{
    U = rUA*(UEqn().H() - fvc::grad(p));
    U.correctBoundaryConditions();
}

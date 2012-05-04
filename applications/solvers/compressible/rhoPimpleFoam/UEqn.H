// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  - fvm::Sp(fvc::ddt(rho) + fvc::div(phi), U)
  + turbulence->divDevRhoReff(U)
);

UEqn().relax();

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));
    K = 0.5*magSqr(U);
}

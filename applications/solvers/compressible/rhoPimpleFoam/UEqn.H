// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
);

UEqn().relax();

sources.constrain(UEqn());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p) + sources(rho, U));
    K = 0.5*magSqr(U);
}

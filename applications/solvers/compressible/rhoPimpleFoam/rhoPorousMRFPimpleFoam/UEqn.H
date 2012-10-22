// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    //pZones.ddt(rho, U)
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
);

UEqn().relax();

mrfZones.addCoriolis(rho, UEqn());

pZones.addResistance(UEqn());

sources.constrain(UEqn());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p) + sources(rho, U));
}

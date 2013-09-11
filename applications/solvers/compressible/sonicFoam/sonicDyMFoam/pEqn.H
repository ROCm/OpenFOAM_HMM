rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

surfaceScalarField phid
(
    "phid",
    fvc::interpolate(psi)
   *(
        (
            (fvc::interpolate(rho*HbyA) & mesh.Sf())
          //***HGW + rhorAUf*fvc::ddtCorr(rho, U, phi)
        )/fvc::interpolate(rho)
      - fvc::meshPhi(rho, U)
    )
);

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
      + fvm::div(phid, p)
      - fvm::laplacian(rhorAUf, p)
    );

    pEqn.solve();

    phi = pEqn.flux();
}

#include "compressibleContinuityErrs.H"

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();

rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField Dp("Dp", fvc::interpolate(rho*rAU));

volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

surfaceScalarField phid
(
    "phid",
    fvc::interpolate(psi)*
    (
        (mesh.Sf() & fvc::interpolate(rho*HbyA))
      + Dp*fvc::ddtCorr(rho, U, phi)
    )/fvc::interpolate(rho)
);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
      + fvm::div(phid, p)
      - fvm::laplacian(Dp, p)
    );

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = pEqn.flux();
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
K = 0.5*magSqr(U);

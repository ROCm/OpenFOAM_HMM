surfaceScalarField Dp("Dp", fvc::interpolate(rAU));

volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn().H();

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(HbyA) & mesh.Sf())
  + Dp*fvc::ddtCorr(U, Uf)
);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(Dp, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

{
    Uf = fvc::interpolate(U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Uf += mesh.Sf()*(phi - (mesh.Sf() & Uf))/sqr(mesh.magSf());
}

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

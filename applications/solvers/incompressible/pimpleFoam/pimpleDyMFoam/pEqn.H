U = rAU*UEqn().H();

if (pimple.nCorr() <= 1)
{
    UEqn.clear();
}

phi = (fvc::interpolate(U) & mesh.Sf());

if (p.needReference())
{
    fvc::makeRelative(phi, U);
    adjustPhi(phi, U, p);
    fvc::makeAbsolute(phi, U);
}

for (int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve
    (
        mesh.solver(p.select(pimple.finalInnerIter(corr, nonOrth)))
    );

    if (nonOrth == pimple.nNonOrthCorr())
    {
        phi -= pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();

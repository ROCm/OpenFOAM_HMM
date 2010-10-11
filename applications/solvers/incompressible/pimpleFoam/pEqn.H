U = rAU*UEqn().H();

if (nCorr <= 1)
{
    UEqn.clear();
}

phi = (fvc::interpolate(U) & mesh.Sf())
    + fvc::ddtPhiCorr(rAU, U, phi);

adjustPhi(phi, U, p);

// Non-orthogonal pressure corrector loop
for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve
    (
        mesh.solver
        (
            p.select
            (
                (
                    finalIter
                 && corr == nCorr-1
                 && nonOrth == nNonOrthCorr
                )
            )
        )
    );

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();

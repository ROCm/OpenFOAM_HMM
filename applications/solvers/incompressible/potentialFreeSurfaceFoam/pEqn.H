volScalarField rAU(1.0/UEqn().A());
surfaceScalarField rAUf(rAU.name() + 'f', fvc::interpolate(rAU));

volVectorField HbyA("HbyA", U);
HbyA = rAU*(UEqn() == sources(U))().H();

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(HbyA) & mesh.Sf())
  + fvc::ddtPhiCorr(rAU, U, phi)
);

adjustPhi(phiHbyA, U, p_gh);


// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_ghEqn
    (
        fvm::laplacian(rAUf, p_gh) == fvc::div(phiHbyA)
    );

    p_ghEqn.setReference(p_ghRefCell, p_ghRefValue);

    p_ghEqn.solve(mesh.solver(p_gh.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - p_ghEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p_gh.relax();

p = p_gh + (g & (mesh.C() + zeta - refLevel));

U = HbyA - rAU*fvc::grad(p_gh);
U.correctBoundaryConditions();
sources.correct(U);

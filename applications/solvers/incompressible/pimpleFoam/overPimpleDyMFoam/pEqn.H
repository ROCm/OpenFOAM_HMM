volScalarField rAU(1.0/UEqn.A());

// Option 1: interpolate rAU, do not block out rAU on blocked cells
//mesh.interpolate(rAU, false);
//surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

// Option 2: do not interpolate rAU but block out rAU
//surfaceScalarField rAUf("rAUf", fvc::interpolate(blockedCells*rAU));


// Option 3: do not interpolate rAU but zero out rAUf on faces on holes
// But what about:
//
//   H
// H I C C C C
//   H
//

surfaceScalarField faceMask(localMin<scalar>(mesh).interpolate(cellMask));
surfaceScalarField rAUf("rAUf", faceMask*fvc::interpolate(rAU));

volVectorField HbyA("HbyA", U);
HbyA = constrainHbyA(cellMask*rAU*UEqn.H(), U, p);

//mesh.interpolate(HbyA);
if (massFluxInterpolation)
{
    #include "interpolatedFaces.H"
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

if (ddtCorr)
{
    phiHbyA += rAUf*fvc::ddtCorr(U, Uf);
}

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

if (adjustFringe)
{
    fvc::makeRelative(phiHbyA, U);
    oversetAdjustPhi(phiHbyA, U);
    fvc::makeAbsolute(phiHbyA, U);
}

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
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
volVectorField gradP(fvc::grad(p));
//mesh.interpolate(gradP);

// Option 1: leave velocity intact on blocked out cells
//U = HbyA - rAU*gradP;

// Option 2: zero out velocity on blocked out cells
U = (HbyA - rAU*cellMask*gradP);
U.correctBoundaryConditions();

fvOptions.correct(U);

{
    Uf = fvc::interpolate(U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Uf += n*(phi/mesh.magSf() - (n & Uf));
}

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

// Option 1: interpolate rAU, do not block out rAU on blocked cells
volScalarField rAU("rAU", 1.0/UEqn.A());
mesh.interpolate(rAU);

// Option 2: do not interpolate rAU but block out rAU
//surfaceScalarField rAUf("rAUf", fvc::interpolate(blockedCells*rAU));


// Option 3: do not interpolate rAU but zero out rAUf on faces on holes
// But what about:
//
//   H
// H I C C C C
//   H
//
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField H("H", UEqn.H());

volVectorField HbyA("HbyA", U);
HbyA = constrainHbyA(rAU*H, U, p);

if (massFluxInterpolation)
{
    #include "interpolatedFaces.H"
}

if (runTime.outputTime())
{
    H.write();
    rAU.write();
    HbyA.write();
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

phiHbyA = fvc::flux(HbyA);

if (ddtCorr)
{
    surfaceScalarField faceMaskOld
    (
        localMin<scalar>(mesh).interpolate(cellMask.oldTime())
    );
    phiHbyA += rAUf*faceMaskOld*fvc::ddtCorr(U, Uf);
}

MRF.makeRelative(phiHbyA);

// WIP
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
        // option 2:
        // rAUf*fvc::snGrad(p)*mesh.magSf();
    }
}


#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();
volVectorField gradP(fvc::grad(p));

// Option 2: zero out velocity on blocked out cells
//U = HbyA - rAU*cellMask*gradP;
// Option 3: zero out velocity on blocked out cells
// This is needed for the scalar Eq (k,epsilon, etc)
// which can use U as source term
U = cellMask*(HbyA - rAU*gradP);
U.correctBoundaryConditions();

fvOptions.correct(U);

{
    Uf = fvc::interpolate(U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Uf += n*(phi/mesh.magSf() - (n & Uf));
}

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

surfaceScalarField faceMask
(
    localMin<scalar>(mesh).interpolate(cellMask)
);
phi *= faceMask;

volScalarField rUA = 1.0/UEqn().A();
U = rUA*UEqn().H();
UEqn.clear();

phi = fvc::interpolate(U) & mesh.Sf();
adjustPhi(phi, U, pd);
phi += fvc::interpolate(beta*gh*rUA)*fvc::snGrad(T)*mesh.magSf();

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pdEqn
    (
        fvm::laplacian(rUA, pd) == fvc::div(phi)
    );

    pdEqn.setReference(pdRefCell, pdRefValue);

    // retain the residual from the first iteration
    if (nonOrth == 0)
    {
        eqnResidual = pdEqn.solve().initialResidual();
        maxResidual = max(eqnResidual, maxResidual);
    }
    else
    {
        pdEqn.solve();
    }

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pdEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
pd.relax();

U -= rUA*(fvc::grad(pd) - beta*gh*fvc::grad(T));
U.correctBoundaryConditions();

rho = thermo->rho();

volScalarField rUA = 1.0/UEqn().A();
U = rUA*UEqn().H();
UEqn.clear();

bool closedVolume = false;

if (transonic)
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(thermo->psi())*(fvc::interpolate(U) & mesh.Sf())
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn0
        (
            fvm::div(phid, p)
          - fvm::laplacian(rho*rUA, p)
        );
        fvScalarMatrix pEqn = pEqn0;
        pEqn.relax(mesh.relaxationFactor("pEqn"));

        pEqn.setReference(pRefCell, pRefValue);

        // retain the residual from the first iteration
        if (nonOrth == 0)
        {
            eqnResidual = pEqn.solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn0.flux();
        }
    }
}
else
{
    phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
    //phi = fvc::interpolate(rho*U) & mesh.Sf();
    closedVolume = adjustPhi(phi, U, p);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rho*rUA, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        // retain the residual from the first iteration
        if (nonOrth == 0)
        {
            eqnResidual = pEqn.solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }
}


#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

rho = thermo->rho();
rho.relax();
Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();

bound(p, pMin);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(thermo->psi()*p))
        /fvc::domainIntegrate(thermo->psi());
}

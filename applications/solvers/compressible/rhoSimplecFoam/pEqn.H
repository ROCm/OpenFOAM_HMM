volScalarField p0 = p;

volScalarField AU = UEqn().A();
volScalarField AtU = AU - UEqn().H1();
U = UEqn().H()/AU;
UEqn.clear();

bool closedVolume = false;

if (transonic)
{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(psi*U) & mesh.Sf()
        );

        surfaceScalarField phic
        (
            "phic",
            fvc::interpolate(rho/AtU - rho/AU)*fvc::snGrad(p)*mesh.magSf()
          + phid*(fvc::interpolate(p) - fvc::interpolate(p, "UD"))
        );

        refCast<mixedFvPatchScalarField>(p.boundaryField()[1]).refValue()
            = p.boundaryField()[1];

        fvScalarMatrix pEqn
        (
            fvm::div(phid, p)
          + fvc::div(phic)
          - fvm::Sp(fvc::div(phid), p)
          + fvc::div(phid)*p
          - fvm::laplacian(rho/AtU, p)
        );

        pEqn.setReference(pRefCell, pRefValue);

        // Retain the residual from the first iteration
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
            phi == phic + pEqn.flux();
        }
    }
}
else
{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        phi = fvc::interpolate(rho*U) & mesh.Sf();
        closedVolume = adjustPhi(phi, U, p);
        phi += fvc::interpolate(rho/AtU - rho/AU)*fvc::snGrad(p)*mesh.magSf();

        fvScalarMatrix pEqn
        (
            fvc::div(phi)
        //- fvm::laplacian(rho/AU, p)
          - fvm::laplacian(rho/AtU, p)
        );

        pEqn.setReference(pRefCell, pRefValue);

        // Retain the residual from the first iteration
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
            phi += pEqn.flux();
        }
    }
}

// The incompressibe for of the continuity error check is appropriate for
// steady-state compressible also.
#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U -= (fvc::grad(p0)*(1.0/AU - 1.0/AtU) + fvc::grad(p)/AtU);
//U -= fvc::grad(p)/AU;

U.correctBoundaryConditions();

bound(p, pMin);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}

rho = thermo.rho();

if (!transonic)
{
    rho.relax();
}

Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

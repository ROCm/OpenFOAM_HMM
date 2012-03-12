rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

volScalarField rAU(1.0/UEqn().A());
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn().H();

UEqn.clear();

bool closedVolume = false;

if (simple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)*(fvc::interpolate(HbyA) & mesh.Sf())
    );

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::div(phid, p)
          - fvm::laplacian(rho*rAU, p)
        );

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax(mesh.equationRelaxationFactor("pEqn"));

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*(fvc::interpolate(HbyA) & mesh.Sf())
    );

    closedVolume = adjustPhi(phiHbyA, U, p);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rho*rAU, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }
}


#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}

rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

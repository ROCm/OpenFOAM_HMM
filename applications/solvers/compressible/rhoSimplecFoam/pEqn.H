rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

volScalarField p0(p);

volScalarField AU(UEqn().A());
volScalarField AtU(AU - UEqn().H1());

volVectorField HbyA("HbyA", U);
HbyA = UEqn().H()/AU;

UEqn.clear();

bool closedVolume = false;

if (simple.transonic())
{
    while (simple.correctNonOrthogonal())
    {
        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(psi*HbyA) & mesh.Sf()
        );

        surfaceScalarField phic
        (
            "phic",
            fvc::interpolate(rho/AtU - rho/AU)*fvc::snGrad(p)*mesh.magSf()
          + phid*(fvc::interpolate(p) - fvc::interpolate(p, "UD"))
        );

        //refCast<mixedFvPatchScalarField>(p.boundaryField()[1]).refValue()
        //    = p.boundaryField()[1];

        fvScalarMatrix pEqn
        (
            fvm::div(phid, p)
          + fvc::div(phic)
          - fvm::Sp(fvc::div(phid), p)
          + fvc::div(phid)*p
          - fvm::laplacian(rho/AtU, p)
        );
        //pEqn.relax();

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi == phic + pEqn.flux();
        }
    }
}
else
{
    while (simple.correctNonOrthogonal())
    {
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::interpolate(rho*HbyA) & mesh.Sf()
        );

        closedVolume = adjustPhi(phi, U, p);
        phi += fvc::interpolate(rho/AtU - rho/AU)*fvc::snGrad(p)*mesh.magSf();

        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
        //- fvm::laplacian(rho/AU, p)
          - fvm::laplacian(rho/AtU, p)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

// The incompressibe for of the continuity error check is appropriate for
// steady-state compressible also.
#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - (fvc::grad(p0)*(1.0/AU - 1.0/AtU) + fvc::grad(p)/AtU);
//U = HbyA - fvc::grad(p)/AU;

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

if (!simple.transonic())
{
    rho.relax();
}

Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

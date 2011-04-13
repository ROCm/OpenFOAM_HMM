rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("(rho*(1|A(U)))", fvc::interpolate(rho*rAU));
U = rAU*UEqn.H();

surfaceScalarField phiU
(
    fvc::interpolate(rho)
   *(
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    )
);

phi = phiU - rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf();

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    surfaceScalarField rhorAUf(fvc::interpolate(rho*rAU));

    fvScalarMatrix p_rghEqn
    (
        fvm::ddt(psi, p_rgh) + fvc::ddt(psi, rho)*gh
      + fvc::div(phi)
      - fvm::laplacian(rhorAUf, p_rgh)
    );

    p_rghEqn.solve
    (
        mesh.solver
        (
            p_rgh.select
            (
                pimple.finalIter()
             && corr == nCorr-1
             && nonOrth == nNonOrthCorr
            )
        )
    );

    if (nonOrth == nNonOrthCorr)
    {
        phi += p_rghEqn.flux();
    }
}

p = p_rgh + rho*gh;

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U += rAU*fvc::reconstruct((phi - phiU)/rhorAUf);
U.correctBoundaryConditions();

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

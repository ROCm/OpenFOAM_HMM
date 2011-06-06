rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf(rAU.name() + 'f', fvc::interpolate(rho*rAU));
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

for (int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++)
{
    fvScalarMatrix p_rghEqn
    (
        fvc::ddt(psi, rho)*gh
      + fvc::div(phi)
      + fvm::ddt(psi, p_rgh)
      - fvm::laplacian(rhorAUf, p_rgh)
     ==
        parcels.Srho()
      + surfaceFilm.Srho()
    );

    p_rghEqn.solve
    (
        mesh.solver(p_rgh.select(pimple.finalInnerIter(corr, nonOrth)))
    );

    if (nonOrth == pimple.nNonOrthCorr())
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

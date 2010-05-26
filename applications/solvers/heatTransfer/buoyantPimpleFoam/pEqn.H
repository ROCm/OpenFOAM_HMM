{
    bool closedVolume = p.needReference();

    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p;

    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));

    U = rUA*UEqn.H();

    phi = fvc::interpolate(rho)*
    (
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

    surfaceScalarField buoyancyPhi =
        rhorUAf*fvc::interpolate(rho)*(g & mesh.Sf());
    phi += buoyancyPhi;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phi)
          - fvm::laplacian(rhorUAf, p)
        );

        pEqn.solve
        (
            mesh.solver
            (
                p.select
                (
                    (
                        finalIter
                     && corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                )
            )
        );

        if (nonOrth == nNonOrthCorr)
        {
            // Calculate the conservative fluxes
            phi += pEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U += rUA*fvc::reconstruct((buoyancyPhi + pEqn.flux())/rhorUAf);
            U.correctBoundaryConditions();
        }
    }

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p +=
            (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
        thermo.rho() = psi*p;
        rho += (initialMass - fvc::domainIntegrate(rho))/totalVolume;
    }
}

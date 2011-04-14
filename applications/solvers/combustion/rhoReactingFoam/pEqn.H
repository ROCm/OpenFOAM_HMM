{
    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p;

    volScalarField rAU(1.0/UEqn.A());
    U = rAU*UEqn.H();

    if (pimple.transonic())
    {
        surfaceScalarField phiv
        (
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        );

        phi = fvc::interpolate(rho)*phiv;

        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(thermo.psi())*phiv
        );

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + fvc::div(phi)
          + correction(fvm::ddt(psi, p) + fvm::div(phid, p))
        );

        for (int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                pDDtEqn
              - fvm::laplacian(rho*rAU, p)
            );

            pEqn.solve
            (
                mesh.solver(p.select(pimple.finalInnerIter(corr, nonOrth)))
            );

            if (nonOrth == pimple.nNonOrthCorr())
            {
                phi += pEqn.flux();
            }
        }
    }
    else
    {
        phi =
            fvc::interpolate(rho)
           *(
                (fvc::interpolate(U) & mesh.Sf())
              + fvc::ddtPhiCorr(rAU, rho, U, phi)
            );

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phi)
        );

        for (int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                pDDtEqn
              - fvm::laplacian(rho*rAU, p)
            );

            pEqn.solve
            (
                mesh.solver(p.select(pimple.finalInnerIter(corr, nonOrth)))
            );

            if (nonOrth == pimple.nNonOrthCorr())
            {
                phi += pEqn.flux();
            }
        }
    }

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();

    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
}

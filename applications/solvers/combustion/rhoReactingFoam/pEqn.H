{
    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p;

    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    if (pimple.transonic())
    {
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            (fvc::interpolate(HbyA) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        );

        surfaceScalarField phid("phid", fvc::interpolate(thermo.psi())*phiHbyA);

        phiHbyA *= fvc::interpolate(rho);

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + fvc::div(phiHbyA)
          + correction(fvm::ddt(psi, p) + fvm::div(phid, p))
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                pDDtEqn
              - fvm::laplacian(rho*rAU, p)
            );

            pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }
    else
    {
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::interpolate(rho)
           *(
                (fvc::interpolate(HbyA) & mesh.Sf())
              + fvc::ddtPhiCorr(rAU, rho, U, phi)
            )
        );

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyA)
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                pDDtEqn
              - fvm::laplacian(rho*rAU, p)
            );

            pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    K = 0.5*magSqr(U);

    dpdt = fvc::ddt(p);
}

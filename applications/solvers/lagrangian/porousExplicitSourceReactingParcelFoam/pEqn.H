{
    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p;

    volScalarField rAU(1.0/UEqn.A());
    U = rAU*UEqn.H();

    if (pZones.size() > 0)
    {
        // ddtPhiCorr not well defined for cases with porosity
        phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
    }
    else
    {
        phi =
            fvc::interpolate(rho)
           *(
                (fvc::interpolate(U) & mesh.Sf())
              + fvc::ddtPhiCorr(rAU, rho, U, phi)
            );
    }

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p))
      + fvc::div(phi)
     ==
        parcels.Srho()
      + massSource.SuTot()
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

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();
}

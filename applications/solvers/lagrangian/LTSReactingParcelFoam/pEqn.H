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

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            pDDtEqn
          - fvm::laplacian(rho*rAU, p)
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
            phi += pEqn.flux();
        }
    }

    Info<< "p min/max = " << min(p).value() << ", " << max(p).value() << endl;

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    #include "rhoEqn.H" // NOTE: flux and time scales now inconsistent
    #include "compressibleContinuityErrs.H"

    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();

    rho = thermo.rho();
    rho = max(rho, rhoMin);
    rho = min(rho, rhoMax);

    #include "setPressureWork.H"
}

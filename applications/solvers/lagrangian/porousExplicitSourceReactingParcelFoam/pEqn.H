{
    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p;

    if (pressureImplicitPorosity)
    {
        U = trTU()&UEqn.H();
    }
    else
    {
        U = trAU()*UEqn.H();
    }

    if (transonic)
    {
        surfaceScalarField phiv = fvc::interpolate(U) & mesh.Sf();

        phi = fvc::interpolate(rho)*phiv;

        surfaceScalarField phid("phid", fvc::interpolate(psi)*phiv);

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            tmp<fvScalarMatrix> lapTerm;

            if (pressureImplicitPorosity)
            {
                lapTerm = fvm::laplacian(rho*trTU(), p);
            }
            else
            {
                lapTerm = fvm::laplacian(rho*trAU(), p);
            }

            fvScalarMatrix pEqn
            (
                fvc::ddt(rho) + fvc::div(phi)
              + correction(psi*fvm::ddt(p) + fvm::div(phid, p))
              - lapTerm()
            ==
                parcels.Srho()
              + pointMassSources.Su()
            );

            if
            (
                oCorr == nOuterCorr-1
             && corr == nCorr-1
             && nonOrth == nNonOrthCorr
            )
            {
                pEqn.solve(mesh.solver("pFinal"));
            }
            else
            {
                pEqn.solve();
            }

            if (nonOrth == nNonOrthCorr)
            {
                phi == pEqn.flux();
            }
        }
    }
    else
    {
        phi =
            fvc::interpolate(rho)
           *(
                (fvc::interpolate(U) & mesh.Sf())
//            + fvc::ddtPhiCorr(trAU(), rho, U, phi)
            );

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            tmp<fvScalarMatrix> lapTerm;

            if (pressureImplicitPorosity)
            {
                lapTerm = fvm::laplacian(rho*trTU(), p);
            }
            else
            {
                lapTerm = fvm::laplacian(rho*trAU(), p);
            }

            fvScalarMatrix pEqn
            (
                fvc::ddt(rho) + psi*correction(fvm::ddt(p))
              + fvc::div(phi)
              - lapTerm()
            ==
                parcels.Srho()
              + pointMassSources.Su()
            );

            if
            (
                oCorr == nOuterCorr-1
             && corr == nCorr-1
             && nonOrth == nNonOrthCorr
            )
            {
                pEqn.solve(mesh.solver("pFinal"));
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

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    if (pressureImplicitPorosity)
    {
        U -= trTU()&fvc::grad(p);
    }
    else
    {
        U -= trAU()*fvc::grad(p);
    }
    U.correctBoundaryConditions();

    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
}

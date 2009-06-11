rho = thermo->rho();

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
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(thermo->psi())
       *(
            (fvc::interpolate(U) & mesh.Sf())
//          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
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
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - lapTerm()
         ==
            reactingParcels.Srho()
          + pointMassSources.Su()
        );

        pEqn.solve();

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
          + fvc::ddtPhiCorr(trAU(), rho, U, phi)
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
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - lapTerm()
         ==
            reactingParcels.Srho()
          + pointMassSources.Su()
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }
}

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


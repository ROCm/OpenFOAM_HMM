{
    bool closedVolume = false;

    rho = thermo.rho();

    volScalarField rUA = 1.0/UEqn().A();
    U = rUA*UEqn().H();

    phi =
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
      - fvc::interpolate(rho*rUA*gh)*fvc::snGrad(rho)*mesh.magSf();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvc::ddt(rho)
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, pd)
        );

        if (corr == nCorr-1 && nonOrth == nNonOrthCorr)
        {
            pdEqn.solve(mesh.solver(pd.name() + "Final"));
        }
        else
        {
            pdEqn.solve(mesh.solver(pd.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi += pdEqn.flux();
        }
    }

    // Update pressure field (including bc)
    thermo.p() == pd + rho*gh + pRef;
    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    // Solve continuity
#   include "rhoEqn.H"

    // Update continuity errors
#   include "compressibleContinuityErrors.H"

    // Correct velocity field
    U -= rUA*(fvc::grad(pd) + fvc::grad(rho)*gh);
    U.correctBoundaryConditions();

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        thermo.p() +=
            (
                dimensionedScalar("massIni", dimMass, initialMassf[i])
              - fvc::domainIntegrate(thermo.psi()*thermo.p())
            )/fvc::domainIntegrate(thermo.psi());
        pd == thermo.p() - (rho*gh + pRef);
        rho = thermo.rho();
    }

    // Update thermal conductivity
    K = thermo.Cp()*turb.alphaEff();
}

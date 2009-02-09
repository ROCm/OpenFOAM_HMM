{
    bool closedVolume = pd.needReference();

    rho = thermo.rho();

    volScalarField rUA = 1.0/UEqn().A();
    surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));

    U = rUA*UEqn().H();

    surfaceScalarField phiU
    (
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
    );

    phi = phiU - ghf*fvc::snGrad(rho)*rhorUAf*mesh.magSf();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::ddt(psi, pd)
          + fvc::ddt(psi)*pRef
          + fvc::ddt(psi, rho)*gh
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

    // Correct velocity field
    U += rUA*fvc::reconstruct((phi - phiU)/rhorUAf);
    U.correctBoundaryConditions();

    // Update pressure field (including bc)
    p == pd + rho*gh + pRef;
    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    // Solve continuity
#   include "rhoEqn.H"

    // Update continuity errors
#   include "compressibleContinuityErrors.H"

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (massIni - fvc::domainIntegrate(psi*p))/fvc::domainIntegrate(psi);
        rho = thermo.rho();
    }

    // Update thermal conductivity
    K = thermoFluid[i].Cp()*turb.alphaEff();

    // Update pd (including bc)
    pd == p - (rho*gh + pRef);
}

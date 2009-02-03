{
    bool closedVolume = pd.needReference();

    rho = thermo->rho();

    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rhorUAf = fvc::interpolate(rho*rUA);

    U = rUA*UEqn.H();

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
          - fvm::laplacian(rhorUAf, pd)
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

    U += rUA*fvc::reconstruct((phi - phiU)/rhorUAf);
    U.correctBoundaryConditions();

    p == pd + rho*gh + pRef;
    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(thermo->psi()*p))
            /fvc::domainIntegrate(thermo->psi());
        rho = thermo->rho();
    }

    pd == p - (rho*gh + pRef);
}

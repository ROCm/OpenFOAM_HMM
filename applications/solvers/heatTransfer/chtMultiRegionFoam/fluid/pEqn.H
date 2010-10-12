{
    bool closedVolume = p_rgh.needReference();

    rho = thermo.rho();

    volScalarField rAU = 1.0/UEqn().A();
    surfaceScalarField rhorAUf("(rho*(1|A(U)))", fvc::interpolate(rho*rAU));

    U = rAU*UEqn().H();

    surfaceScalarField phiU
    (
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        )
    );

    phi = phiU - rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::ddt(psi, p_rgh) + fvc::ddt(psi, rho)*gh
          + fvc::div(phi)
          - fvm::laplacian(rhorAUf, p_rgh)
        );

        p_rghEqn.solve
        (
            mesh.solver
            (
                p_rgh.select
                (
                    (
                        oCorr == nOuterCorr-1
                     && corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                )
            )
        );

        if (nonOrth == nNonOrthCorr)
        {
            phi += p_rghEqn.flux();
        }
    }

    // Correct velocity field
    U += rAU*fvc::reconstruct((phi - phiU)/rhorAUf);
    U.correctBoundaryConditions();

    p = p_rgh + rho*gh;

    // Update pressure substantive derivative
    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    // Solve continuity
    #include "rhoEqn.H"

    // Update continuity errors
    #include "compressibleContinuityErrors.H"

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
        rho = thermo.rho();
        p_rgh = p - rho*gh;
    }

    // Update thermal conductivity
    K = thermoFluid[i].Cp()*turb.alphaEff();
}

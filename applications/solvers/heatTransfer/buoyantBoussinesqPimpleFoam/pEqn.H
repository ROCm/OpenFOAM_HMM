{
    volScalarField rUA("rUA", 1.0/UEqn.A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U = rUA*UEqn.H();

    phi = (fvc::interpolate(U) & mesh.Sf())
        + fvc::ddtPhiCorr(rUA, U, phi);

    surfaceScalarField buoyancyPhi =
        rUAf*fvc::interpolate(rhok)*(g & mesh.Sf());
    phi += buoyancyPhi;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUAf, p) == fvc::div(phi)
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
            // Calculate the conservative fluxes
            phi -= pEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U += rUA*fvc::reconstruct((buoyancyPhi - pEqn.flux())/rUAf);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"
}

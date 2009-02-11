{
    volScalarField rUA("rUA", 1.0/UEqn().A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U = rUA*UEqn().H();
    UEqn.clear();

    phi =
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, U, phi)
      + betaghf*fvc::snGrad(T)*rUAf*mesh.magSf();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd) == fvc::div(phi)
        );

        // retain the residual from the first iteration
        if (nonOrth == 0)
        {
            pdEqn.solve(mesh.solver(pd.name() + "Final"));
        }
        else
        {
            pdEqn.solve(mesh.solver(pd.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            // Calculate the conservative fluxes
            phi -= pdEqn.flux();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U -=
                rUA
               *fvc::reconstruct((pdEqn.flux() - betaghf*fvc::snGrad(T))/rUAf);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"
}

{
    volScalarField rUA("rUA", 1.0/UEqn().A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U = rUA*UEqn().H();
    UEqn.clear();

    phi = fvc::interpolate(U) & mesh.Sf();
    adjustPhi(phi, U, pd);
    surfaceScalarField buoyancyPhi = -betaghf*fvc::snGrad(T)*rUAf*mesh.magSf();
    phi -= buoyancyPhi;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd) == fvc::div(phi)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

        // retain the residual from the first iteration
        if (nonOrth == 0)
        {
            eqnResidual = pdEqn.solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            pdEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            // Calculate the conservative fluxes
            phi -= pdEqn.flux();

            // Explicitly relax pressure for momentum corrector
            pd.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U -= rUA*fvc::reconstruct((buoyancyPhi + pdEqn.flux())/rUAf);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"
}

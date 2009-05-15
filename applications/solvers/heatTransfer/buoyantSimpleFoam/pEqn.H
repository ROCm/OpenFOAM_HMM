{
    volScalarField rUA = 1.0/UEqn().A();
    surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));

    U = rUA*UEqn().H();
    UEqn.clear();

    phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
    bool closedVolume = adjustPhi(phi, U, p);
    surfaceScalarField buoyancyPhi = ghf*fvc::snGrad(rho)*rhorUAf*mesh.magSf();
    phi -= buoyancyPhi;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rhorUAf, pd) == fvc::div(phi)
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
            U -= rUA*fvc::reconstruct((buoyancyPhi + pdEqn.flux())/rhorUAf);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"

    p == pd + rho*gh + pRef;

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(thermo->psi()*p))
            /fvc::domainIntegrate(thermo->psi());
    }

    rho = thermo->rho();
    rho.relax();
    Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value()
        << endl;

    pd == p - (rho*gh + pRef);
}

{
    rho = thermo.rho();

    volScalarField rAU(1.0/UEqn().A());
    surfaceScalarField rhorAUf("(rho*(1|A(U)))", fvc::interpolate(rho*rAU));

    U = rAU*UEqn().H();
    UEqn.clear();

    phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
    bool closedVolume = adjustPhi(phi, U, p);

    surfaceScalarField buoyancyPhi
    (
        rhorAUf*fvc::interpolate(rho)*(g & mesh.Sf())
    );
    phi += buoyancyPhi;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rhorAUf, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            // For closed-volume cases adjust the pressure and density levels
            // to obey overall mass continuity
            if (closedVolume)
            {
                p += (initialMass - fvc::domainIntegrate(psi*p))
                    /fvc::domainIntegrate(psi);
            }

            // Calculate the conservative fluxes
            phi -= pEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U += rAU*(rho*g - fvc::grad(p));
            //U += rAU*fvc::reconstruct((buoyancyPhi - pEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"

    rho = thermo.rho();
    rho.relax();
    Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value()
        << endl;
}

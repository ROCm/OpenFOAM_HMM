{
    volScalarField rAU = 1.0/UEqn.A();
    surfaceScalarField rAUf = fvc::interpolate(rAU);

    volVectorField HU = UEqn.H();
    U = rAU*HU;

    surfaceScalarField phiU("phiU", (fvc::interpolate(U) & mesh.Sf()));

    if (ddtPhiCorr)
    {
        phiU += fvc::ddtPhiCorr(rAU, rho, U, phiAbs);
    }

    phi = phiU +
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(gamma)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf();

    if (pd.needReference())
    {
        if (mesh.changing())
        {
            fvc::makeRelative(phi, U);
            adjustPhi(phi, U, pd);
            fvc::makeAbsolute(phi, U);
        }
        else
        {
            adjustPhi(phi, U, pd);
        }
    }

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rAUf, pd) == fvc::div(phi)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

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
            phi -= pdEqn.flux();
        }
    }

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}

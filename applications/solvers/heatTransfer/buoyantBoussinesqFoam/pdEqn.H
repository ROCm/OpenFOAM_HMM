{
    volScalarField rUA("rUA", 1.0/UEqn().A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U = rUA*UEqn().H();
    UEqn.clear();

    surfaceScalarField phiU
    (
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, U, phi)
    );

    phi = phiU + betaghf*fvc::snGrad(T)*rUAf*mesh.magSf();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd) == fvc::div(phi)
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

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();

    #include "continuityErrs.H"
}

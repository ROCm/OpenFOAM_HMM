{
    p.boundaryField().updateCoeffs();

    volScalarField rAUrel(1.0/UrelEqn().A());
    Urel = rAUrel*UrelEqn().H();
    UrelEqn.clear();

    phi = fvc::interpolate(Urel, "interpolate(HbyA)") & mesh.Sf();
    adjustPhi(phi, Urel, p);

    // Non-orthogonal pressure corrector loop
    for (int nonOrth=0; nonOrth<=simple.nNonOrthCorr(); nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUrel, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (nonOrth == simple.nNonOrthCorr())
        {
            phi -= pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    Urel -= rAUrel*fvc::grad(p);
    Urel.correctBoundaryConditions();
}

{
    volVectorField HbyA("HbyA", U);
    HbyA = rAUc*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + alphacf*fvc::interpolate(rho*rAUc)*fvc::ddtCorr(U, phi)
    );

    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        phicForces +
        (
            fvc::interpolate(mixture.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUcf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - MRF.relative(mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUcf.boundaryField())
    );

    while (pimple.correctNonOrthogonal())
    {
        surfaceScalarField Dp("Dp", alphacf*rAUcf);

        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(Dp, p_rgh)
         ==
            fvc::ddt(alphac) + fvc::div(alphacf*phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux()/alphacf;

            p_rgh.relax();

            U =
                HbyA
              + rAUc*fvc::reconstruct((phig - p_rghEqn.flux()/alphacf)/rAUcf);

            U.correctBoundaryConditions();

            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}

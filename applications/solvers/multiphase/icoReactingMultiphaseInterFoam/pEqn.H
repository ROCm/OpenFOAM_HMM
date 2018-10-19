{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        (
            fluid.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    while (pimple.correctNonOrthogonal())
    {

        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        forAllConstIters(fluid.totalPhasePairs(), iter)
        {
            const phasePair& pair = iter()();

            const phaseModel& phase1 = pair.phase1();
            const phaseModel& phase2 = pair.phase2();

            const phasePairKey key12
            (
                phase1.name(),
                phase2.name(),
                true
            );

            // Mass transfer from phase2 to phase1
            tmp<volScalarField> tdmdt12(fluid.dmdt(key12));
            const volScalarField& dmdt12 = tdmdt12();

            const phasePairKey key21
            (
                phase2.name(),
                phase1.name(),
                true
            );

            // Mass transfer from phase1 to phase2
            tmp<volScalarField> tdmdt21(fluid.dmdt(key21));
            const volScalarField&  dmdt21 = tdmdt21();

            const volScalarField dmdtNet(dmdt21 - dmdt12);

            p_rghEqn +=
                dmdtNet*
                (
                  - fluid.coeffs(phase1.name())
                  + fluid.coeffs(phase2.name())
                );
        }

        p_rghEqn.setReference(pRefCell, pRefValue);

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);

            U.correctBoundaryConditions();
            fvOptions.correct(U);
            K = 0.5*magSqr(U);
        }
    }

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

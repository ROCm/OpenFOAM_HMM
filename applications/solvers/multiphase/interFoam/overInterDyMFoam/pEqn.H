{
    rAU = 1.0/UEqn.A();
    //mesh.interpolate(rAU);

    surfaceScalarField faceMask(localMin<scalar>(mesh).interpolate(cellMask));

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField H("H", UEqn.H());

    volVectorField HbyA("HbyA", U);
    //HbyA = rAU*UEqn.H();
    HbyA = constrainHbyA(rAU*H, U, p_rgh);

    if (massFluxInterpolation)
    {
        #include "interpolatedFaces.H"
    }

    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

    if (ddtCorr)
    {
        surfaceScalarField faceMaskOld
        (
            localMin<scalar>(mesh).interpolate(cellMask.oldTime())
        );
        phiHbyA +=
            fvc::interpolate(rho*rAU)*faceMaskOld*fvc::ddtCorr(U, Uf);
    }
    MRF.makeRelative(phiHbyA);

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    if (adjustFringe)
    {
        fvc::makeRelative(phiHbyA, U);
        oversetAdjustPhi(phiHbyA, U);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*faceMask*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    if (adjustFringe)
    {
        fvc::makeRelative(phiHbyA, U);
        oversetAdjustPhi(phiHbyA, U);
        fvc::makeAbsolute(phiHbyA, U);
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U =
                cellMask*
                (
                   HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf)
                );

            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    {
        Uf = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf += n*(phi/mesh.magSf() - (n & Uf));
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    // Zero faces H-I for transport Eq after pEq
    phi *= faceMask;

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p_rgh += dimensionedScalar
        (
            "p_rgh",
            p.dimensions(),
            pRefValue - getRefCellValue(p_rgh, pRefCell)
        );
        p == p_rgh + rho*gh;
    }
}

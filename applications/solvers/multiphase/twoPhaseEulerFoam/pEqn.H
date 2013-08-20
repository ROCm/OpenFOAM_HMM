{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1));
    surfaceScalarField alpha2f(scalar(1) - alpha1f);

    rAU1 = 1.0/U1Eqn.A();
    rAU2 = 1.0/U2Eqn.A();

    surfaceScalarField rAlphaAU1f(fvc::interpolate(alpha1*rAU1));
    surfaceScalarField rAlphaAU2f(fvc::interpolate(alpha2*rAU2));

    // Update the phi BCs from Us before p BCs are updated
    phi.boundaryField() =
        alpha1f.boundaryField()
       *(mesh.Sf().boundaryField() & U1.boundaryField())
      + alpha2f.boundaryField()
       *(mesh.Sf().boundaryField() & U2.boundaryField());

    volVectorField HbyA1
    (
        IOobject::groupName("HbyA", phase1.name()),
        U1
    );
    HbyA1 = rAU1*U1Eqn.H();

    volVectorField HbyA2
    (
        IOobject::groupName("HbyA", phase2.name()),
        U2
    );
    HbyA2 = rAU2*U2Eqn.H();

    mrfZones.absoluteFlux(phi1.oldTime());
    mrfZones.absoluteFlux(phi1);
    mrfZones.absoluteFlux(phi2.oldTime());
    mrfZones.absoluteFlux(phi2);

    // Phase-1 pressure flux (e.g. due to particle-particle pressure)
    surfaceScalarField phiP1
    (
        "phiP1",
        fvc::interpolate((1.0/rho1)*rAU1*turbulence1().pPrime())
       *fvc::snGrad(alpha1)*mesh.magSf()
    );

    // Phase-2 pressure flux (e.g. due to particle-particle pressure)
    surfaceScalarField phiP2
    (
        "phiP2",
        fvc::interpolate((1.0/rho2)*rAU2*turbulence2().pPrime())
       *fvc::snGrad(alpha2)*mesh.magSf()
    );

    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1.name()),
        (fvc::interpolate(HbyA1) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU1, alpha1, U1, phi1)
    );

    surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2.name()),
        (fvc::interpolate(HbyA2) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU2, alpha2, U2, phi2)
    );

    phi = alpha1f*phiHbyA1 + alpha2f*phiHbyA2;
    mrfZones.relativeFlux(phi);

    phiHbyA1 +=
    (
        fvc::interpolate((1.0/rho1)*rAU1*dragCoeff)*phi2
      - phiP1
      + rAlphaAU1f*(g & mesh.Sf())
    );
    mrfZones.relativeFlux(phiHbyA1);

    phiHbyA2 +=
    (
        fvc::interpolate((1.0/rho2)*rAU2*dragCoeff)*phi1
      - phiP2
      + rAlphaAU2f*(g & mesh.Sf())
    );
    mrfZones.relativeFlux(phiHbyA2);

    mrfZones.relativeFlux(phi1.oldTime());
    mrfZones.relativeFlux(phi1);
    mrfZones.relativeFlux(phi2.oldTime());
    mrfZones.relativeFlux(phi2);

    surfaceScalarField phiHbyA("phiHbyA", alpha1f*phiHbyA1 + alpha2f*phiHbyA2);

    HbyA1 += (1.0/rho1)*rAU1*dragCoeff*U2;
    HbyA2 += (1.0/rho2)*rAU2*dragCoeff*U1;

    surfaceScalarField Dp
    (
        "Dp",
        mag
        (
            alpha1f*rAlphaAU1f/fvc::interpolate(rho1)
          + alpha2f*rAlphaAU2f/fvc::interpolate(rho2)
        )
    );

    tmp<fvScalarMatrix> pEqnComp1;
    tmp<fvScalarMatrix> pEqnComp2;

    if (pimple.transonic())
    {
        surfaceScalarField phid1
        (
            IOobject::groupName("phid", phase1.name()),
            fvc::interpolate(psi1)*phi1
        );
        surfaceScalarField phid2
        (
            IOobject::groupName("phid", phase2.name()),
            fvc::interpolate(psi2)*phi2
        );

        pEqnComp1 =
            fvc::ddt(rho1)
          + fvc::div(phi1, rho1) - fvc::Sp(fvc::div(phi1), rho1)
          + correction
            (
                psi1*fvm::ddt(p)
              + fvm::div(phid1, p) - fvm::Sp(fvc::div(phid1), p)
            );
        deleteDemandDrivenData(pEqnComp1().faceFluxCorrectionPtr());
        pEqnComp1().relax();

        pEqnComp2 =
            fvc::ddt(rho2)
          + fvc::div(phi2, rho2) - fvc::Sp(fvc::div(phi2), rho2)
          + correction
            (
                psi2*fvm::ddt(p)
              + fvm::div(phid2, p) - fvm::Sp(fvc::div(phid2), p)
            );
        deleteDemandDrivenData(pEqnComp2().faceFluxCorrectionPtr());
        pEqnComp2().relax();
    }
    else
    {
        pEqnComp1 =
            fvc::ddt(rho1) + psi1*correction(fvm::ddt(p))
          + fvc::div(phi1, rho1) - fvc::Sp(fvc::div(phi1), rho1);

        pEqnComp2 =
            fvc::ddt(rho2) + psi2*correction(fvm::ddt(p))
          + fvc::div(phi2, rho2) - fvc::Sp(fvc::div(phi2), rho2);
    }

    // Cache p prior to solve for density update
    volScalarField p_0(p);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(Dp, p)
        );

        solve
        (
            (
                (alpha1/rho1)*pEqnComp1()
              + (alpha2/rho2)*pEqnComp2()
            )
          + pEqnIncomp,
            mesh.solver(p.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField mSfGradp(pEqnIncomp.flux()/Dp);

            phi1.boundaryField() ==
                (mesh.Sf().boundaryField() & U1.boundaryField());
            phi1 = phiHbyA1 + rAlphaAU1f*mSfGradp/fvc::interpolate(rho1);

            phi2.boundaryField() ==
                (mesh.Sf().boundaryField() & U2.boundaryField());
            phi2 = phiHbyA2 + rAlphaAU2f*mSfGradp/fvc::interpolate(rho2);

            phi = alpha1f*phi1 + alpha2f*phi2;

            dgdt =
            (
                pos(alpha2)*(pEqnComp2 & p)/rho2
              - pos(alpha1)*(pEqnComp1 & p)/rho1
            );

            p.relax();
            mSfGradp = pEqnIncomp.flux()/Dp;

            U1 = HbyA1
              + fvc::reconstruct
                (
                    rAlphaAU1f
                   *(
                        (g & mesh.Sf())
                      + mSfGradp/fvc::interpolate(rho1)
                    )
                  - phiP1
                );
            U1.correctBoundaryConditions();

            U2 = HbyA2
              + fvc::reconstruct
                (
                    rAlphaAU2f
                   *(
                        (g & mesh.Sf())
                      + mSfGradp/fvc::interpolate(rho2)
                    )
                  - phiP2
                );
            U2.correctBoundaryConditions();

            U = fluid.U();
        }
    }

    p = max(p, pMin);

    // Update densities from change in p
    rho1 += psi1*(p - p_0);
    rho2 += psi2*(p - p_0);

    K1 = 0.5*magSqr(U1);
    K2 = 0.5*magSqr(U2);

    if (thermo1.dpdt() || thermo2.dpdt())
    {
        dpdt = fvc::ddt(p);
    }
}

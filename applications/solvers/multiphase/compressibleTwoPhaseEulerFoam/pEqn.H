{
    rho1 = rho10 + psi1*p;
    rho2 = rho20 + psi2*p;

    surfaceScalarField alpha1f(fvc::interpolate(alpha1));
    surfaceScalarField alpha2f(scalar(1) - alpha1f);

    volScalarField rAU1(1.0/U1Eqn.A());
    volScalarField rAU2(1.0/U2Eqn.A());

    surfaceScalarField rAlphaAU1f = fvc::interpolate(alpha1*rAU1);
    surfaceScalarField rAlphaAU2f = fvc::interpolate(alpha2*rAU2);

    volVectorField HbyA1("HbyA1", U1);
    HbyA1 = rAU1*U1Eqn.H();

    volVectorField HbyA2("HbyA2", U2);
    HbyA2 = rAU2*U2Eqn.H();

    surfaceScalarField phiHbyA1
    (
        "phiHbyA1",
        (fvc::interpolate(HbyA1) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU1, alpha1, U1, phi1)
      + fvc::interpolate((1.0/rho1)*rAU1*dragCoeff)*phi2
      + rAlphaAU1f*(g & mesh.Sf())
    );

    surfaceScalarField phiHbyA2
    (
        "phiHbyA2",
        (fvc::interpolate(HbyA2) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU2, alpha2, U2, phi2)
      + fvc::interpolate((1.0/rho2)*rAU2*dragCoeff)*phi1
      + rAlphaAU2f*(g & mesh.Sf())
    );

    surfaceScalarField phiHbyA("phiHbyA", alpha1f*phiHbyA1 + alpha2f*phiHbyA2);

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

    //if (transonic)
    //{
    //}
    //else
    {
        surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi1);
        surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi2);

        pEqnComp1 =
            fvc::ddt(rho1) + psi1*correction(fvm::ddt(p))
          + fvc::div(phid1, p)
          - fvc::Sp(fvc::div(phid1), p);

        pEqnComp2 =
            fvc::ddt(rho2) + psi2*correction(fvm::ddt(p))
          + fvc::div(phid2, p)
          - fvc::Sp(fvc::div(phid2), p);
    }

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
            surfaceScalarField mSfGradp = pEqnIncomp.flux()/Dp;

            phi1.boundaryField() ==
                (fvc::interpolate(U1) & mesh.Sf())().boundaryField();
            phi1 = phiHbyA1 + rAlphaAU1f*mSfGradp/fvc::interpolate(rho1);

            phi2.boundaryField() ==
                (fvc::interpolate(U2) & mesh.Sf())().boundaryField();
            phi2 = phiHbyA2 + rAlphaAU2f*mSfGradp/fvc::interpolate(rho2);

            phi = alpha1f*phi1 + alpha2f*phi2;

            dgdt =
            (
                pos(alpha2)*(pEqnComp2 & p)/rho2
              - pos(alpha1)*(pEqnComp1 & p)/rho1
            );

            p.relax();
            mSfGradp = pEqnIncomp.flux()/Dp;

            volVectorField U10 = U1;

            U1 = HbyA1
               + (1.0/rho1)*rAU1*dragCoeff*U2
               + fvc::reconstruct
                 (
                     rAlphaAU1f*(g & mesh.Sf())
                   + rAlphaAU1f*mSfGradp/fvc::interpolate(rho1)
                 );
            U1.correctBoundaryConditions();

            U2 = HbyA2
               + (1.0/rho2)*rAU2*dragCoeff*U10
               + fvc::reconstruct
                 (
                     rAlphaAU2f*(g & mesh.Sf())
                   + rAlphaAU2f*mSfGradp/fvc::interpolate(rho2)
                 );
            U2.correctBoundaryConditions();

            U = alpha1*U1 + alpha2*U2;
        }
    }

    p = max(p, pMin);

    rho1 = rho10 + psi1*p;
    rho2 = rho20 + psi2*p;

    K1 = 0.5*magSqr(U1);
    K2 = 0.5*magSqr(U1);

    dpdt = fvc::ddt(p);
}

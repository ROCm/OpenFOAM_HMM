{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> tphiAlpha;

    if (MULESCorr)
    {
        mesh.interpolate(alpha1);

        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
        );

        alpha1Eqn.solve();

        tmp<surfaceScalarField> tphiAlphaUD(alpha1Eqn.flux());
        alphaPhi = tphiAlphaUD();

        if (alphaApplyPrevCorr && tphiAlphaCorr0.valid())
        {
            MULES::correct(alpha1, alphaPhi, tphiAlphaCorr0.ref(), 1, 0);

            alphaPhi += tphiAlphaCorr0();
        }

        // Cache the upwind-flux
        tphiAlphaCorr0 = tphiAlphaUD;

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phir(phic*mixture.nHatf());

        tmp<surfaceScalarField> tphiAlphaUn
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tphiAlphaCorr(tphiAlphaUn() - alphaPhi);
            volScalarField alpha10("alpha10", alpha1);

            MULES::correct(alpha1, tphiAlphaUn(), tphiAlphaCorr.ref(), 1, 0);
            mesh.interpolate(alpha1);

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi += tphiAlphaCorr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                alphaPhi += 0.5*tphiAlphaCorr();
            }
        }
        else
        {
            alphaPhi = tphiAlphaUn;
            MULES::explicitSolve(alpha1, phi, alphaPhi, 1, 0);
            mesh.interpolate(alpha1);
        }

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }


    rhoPhi = alphaPhi*(rho1 - rho2) + phi*rho2;

    if (alphaApplyPrevCorr && MULESCorr)
    {
        tphiAlphaCorr0 = alphaPhi - tphiAlphaCorr0;
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}

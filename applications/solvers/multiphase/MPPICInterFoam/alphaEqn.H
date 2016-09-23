{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // Standard face-flux compression coefficient
    surfaceScalarField phic
    (
        mixture.cAlpha()*mag(alphaPhic/mesh.magSf())
    );

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    surfaceScalarField::Boundary& phicBf = phic.boundaryFieldRef();
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
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alphac, alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                alphaPhic,
                upwind<scalar>(mesh, alphaPhic)
            ).fvmDiv(alphaPhic, alpha1)
          - fvm::Sp(fvc::ddt(alphac) + fvc::div(alphaPhic), alpha1)
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;

        tmp<surfaceScalarField> tphiAlphaUD(alpha1Eqn.flux());
        alphaPhi = tphiAlphaUD();

        if (alphaApplyPrevCorr && tphiAlphaCorr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;

            MULES::correct
            (
                alphac,
                alpha1,
                alphaPhi,
                tphiAlphaCorr0.ref(),
                zeroField(), zeroField(),
                1, 0
            );

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
                alphaPhic,
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

            //MULES::correct(alpha1, tphiAlphaUn(), tphiAlphaCorr(), 1, 0);

            MULES::correct
            (
                alphac,
                alpha1,
                tphiAlphaUn(),
                tphiAlphaCorr.ref(),
                zeroField(), zeroField(),
                1, 0
            );

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

            MULES::explicitSolve
            (
                alphac,
                alpha1,
                alphaPhic,
                alphaPhi,
                zeroField(), zeroField(),
                1, 0
            );

        }

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }

    rhoPhi = alphaPhi*(rho1 - rho2) + alphaPhic*rho2;

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

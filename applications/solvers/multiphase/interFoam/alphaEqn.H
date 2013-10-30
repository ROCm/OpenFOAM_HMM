{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir(phic*interface.nHatf());

    if (pimple.firstIter() && MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            #ifdef LTSSOLVE
            fv::localEulerDdtScheme<scalar>(mesh, rDeltaT.name()).fvmDdt(alpha1)
            #else
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
            #endif
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;

        tphiAlpha = alpha1Eqn.flux();
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        tmp<surfaceScalarField> tphiAlpha0
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
            tphiAlpha0() -= tphiAlpha();
            #ifdef LTSSOLVE
            MULES::LTScorrect(alpha1, tphiAlpha0(), 1, 0);
            #else
            MULES::correct(alpha1, tphiAlpha0(), 1, 0);
            #endif
            tphiAlpha() += tphiAlpha0();
        }
        else
        {
            tphiAlpha = tphiAlpha0;

            #ifdef LTSSOLVE
            MULES::explicitLTSSolve(alpha1, phi, tphiAlpha(), 1, 0);
            #else
            MULES::explicitSolve(alpha1, phi, tphiAlpha(), 1, 0);
            #endif
        }

        alpha2 = 1.0 - alpha1;
    }

    rhoPhi = tphiAlpha()*(rho1 - rho2) + phi*rho2;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}

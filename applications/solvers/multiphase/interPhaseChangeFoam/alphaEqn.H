{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir("phir", phic*interface.nHatf());

    Pair<tmp<volScalarField> > vDotAlphal =
        twoPhaseProperties->vDotAlphal();
    const volScalarField& vDotcAlphal = vDotAlphal[0]();
    const volScalarField& vDotvAlphal = vDotAlphal[1]();
    const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);

    tmp<surfaceScalarField> tphiAlpha;

    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(alpha1)
          + fvm::div(phi, alpha1, "UD") - fvm::Sp(divU, alpha1)
         ==
            fvm::Sp(vDotvmcAlphal, alpha1)
          + vDotcAlphal
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;

        tphiAlpha = alpha1Eqn.flux();
    }

    volScalarField alpha10("alpha10", alpha1);

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        tmp<surfaceScalarField> tphiAlphaCorr
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
            tphiAlphaCorr() -= tphiAlpha();

            volScalarField alpha100("alpha100", alpha10);
            alpha10 = alpha1;

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                tphiAlphaCorr(),
                vDotvmcAlphal,
                (
                    divU*(alpha10 - alpha100)
                  - vDotvmcAlphal*alpha10
                )(),
                1,
                0
            );

            tphiAlpha() += tphiAlphaCorr();
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                tphiAlphaCorr(),
                vDotvmcAlphal,
                (divU*alpha1 + vDotcAlphal)(),
                1,
                0
            );

            tphiAlpha = tphiAlphaCorr;
        }

        alpha2 = 1.0 - alpha1;
    }

    rhoPhi = tphiAlpha()*(rho1 - rho2) + phi*rho2;

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}

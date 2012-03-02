{
    word scheme("div(phi,alpha1)");
    word schemer("div(phir,alpha1)");

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phi1 - phi2);

    if (g0.value() > 0.0)
    {
        surfaceScalarField alpha1f(fvc::interpolate(alpha1));
        surfaceScalarField phipp(ppMagf*fvc::snGrad(alpha1)*mesh.magSf());
        phir += phipp;
        phic += fvc::interpolate(alpha1)*phipp;
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        fvScalarMatrix alpha1Eqn
        (
             fvm::ddt(alpha1)
           + fvm::div(phic, alpha1, scheme)
           + fvm::div(-fvc::flux(-phir, alpha2, schemer), alpha1, schemer)
        );

        if (g0.value() > 0.0)
        {
            ppMagf = rAU1f*fvc::interpolate
            (
                (1.0/(rho1*(alpha1 + scalar(0.0001))))
               *g0*min(exp(preAlphaExp*(alpha1 - alphaMax)), expMax)
            );

            alpha1Eqn -= fvm::laplacian
            (
                (fvc::interpolate(alpha1) + scalar(0.0001))*ppMagf,
                alpha1,
                "laplacian(alpha1PpMag,alpha1)"
            );
        }

        alpha1Eqn.relax();
        alpha1Eqn.solve();

        #include "packingLimiter.H"

        alpha2 = scalar(1) - alpha1;

        Info<< "Dispersed phase volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;
    }
}

rho = alpha1*rho1 + alpha2*rho2;

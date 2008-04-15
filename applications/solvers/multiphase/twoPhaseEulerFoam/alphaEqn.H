{
    word scheme("div(phi,alpha)");
    word schemer("div(phir,alpha)");

    surfaceScalarField phic = phi;
    surfaceScalarField phir = phia - phib;

    if (g0.value() > 0.0)
    {
        surfaceScalarField alphaf = fvc::interpolate(alpha);
        surfaceScalarField phipp = ppMagf*fvc::snGrad(alpha)*mesh.magSf();
        phir += phipp;
        phic += fvc::interpolate(alpha)*phipp;
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        fvScalarMatrix alphaEqn
        (
             fvm::ddt(alpha)
           + fvm::div(phic, alpha, scheme)
           + fvm::div(-fvc::flux(-phir, beta, schemer), alpha, schemer)
        );

        if (g0.value() > 0.0)
        {
            surfaceScalarField alphaf = fvc::interpolate(alpha);

            // Correct the particle-particle force magnitude
            ppMagf = g0*rUaAf*min(exp(preAlphaExp*(alphaf - alphaMax)), expMax);

            alphaEqn -= fvm::laplacian(ppMagf, alpha);
        }

        alphaEqn.relax();
        alphaEqn.solve();

#       include "packingLimiter.H"

        beta = scalar(1) - alpha;

        Info<< "Dispersed phase volume fraction = "
            << alpha.weightedAverage(mesh.V()).value()
            << "  Min(alpha) = " << min(alpha).value()
            << "  Max(alpha) = " << max(alpha).value()
            << endl;
    }
}

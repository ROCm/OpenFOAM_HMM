{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix hEqn
        (
           - fvm::laplacian(betav*alpha, h, "laplacian(alpha,h)")
           + sources(rho, h)
        );

        hEqn.relax();

        sources.constrain(hEqn);

        hEqn.solve();

        sources.correct(h);
    }
}

thermo.correct();

Info<< "Min/max T:" << min(thermo.T()) << ' ' << max(thermo.T()) << endl;

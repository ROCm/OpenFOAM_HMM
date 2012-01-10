{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        tmp<fvScalarMatrix> TEqn
        (
          - fvm::laplacian(kappa, T, "laplacian(K,T)")
          + sources(rho, T)
        );

        TEqn().relax();
        TEqn().solve();
    }

    Info<< "Min/max T:" << min(T).value() << ' ' << max(T).value() << endl;
}

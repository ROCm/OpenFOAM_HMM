{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        solve
        (
            fvm::ddt(rho*cp, T) - fvm::laplacian(K, T)
        );
    }

    Info<< "Min/max T:" << min(T) << ' ' << max(T) << endl;
}

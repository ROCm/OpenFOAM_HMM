if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        tmp<fvScalarMatrix> TEqn
        (
            fvm::ddt(betav*rho*cp, T)
          - fvm::laplacian(betav*kappa, T, "laplacian(K,T)")
          + sources(rho, T)
        );

        TEqn().relax();
        TEqn().solve(mesh.solver(T.select(finalIter)));
    }

    Info<< "Min/max T:" << min(T).value() << ' ' << max(T).value() << endl;
}

thermo.correct();

if (finalIter)
{
    mesh.data::remove("finalIteration");
}

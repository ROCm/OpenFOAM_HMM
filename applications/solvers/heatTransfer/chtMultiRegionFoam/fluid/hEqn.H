{
    tmp<fvScalarMatrix> hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(phi, h)
      - fvm::laplacian(turb.alphaEff(), h)
     ==
        DpDt
    );
    hEqn().relax();
    hEqn().solve();

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()) << ' ' << max(thermo.T())
        << endl;
}

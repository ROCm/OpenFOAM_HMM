{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(phi, h)
      - fvm::laplacian(turb.alphaEff(), h)
     ==
        DpDt
      + rad.Sh(thermo)
    );

    hEqn.relax();
    hEqn.solve(mesh.solver(h.select(finalIter)));

    thermo.correct();

    rad.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}

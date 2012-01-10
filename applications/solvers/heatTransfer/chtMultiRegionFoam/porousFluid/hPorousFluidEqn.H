{
    fvScalarMatrix hPorousEqn
    (
        fvm::ddt(porousRho, porousH)
      + fvm::div(porousPhi, porousH)
      - fvm::laplacian(turbPorous.alphaEff(), porousH)
     ==
        porousdpdt
      - (fvc::ddt(porousRho, porousK) + fvc::div(porousPhi, porousK))
      + porousSources(porousRho, porousH)
    );

    hPorousEqn.relax();
    hPorousEqn.solve();

    porousThermo.correct();

    Info<< "Min/max in the porous T:"
        << min(porousThermo.T()).value() << ' '
        << max(porousThermo.T()).value() << endl;
}

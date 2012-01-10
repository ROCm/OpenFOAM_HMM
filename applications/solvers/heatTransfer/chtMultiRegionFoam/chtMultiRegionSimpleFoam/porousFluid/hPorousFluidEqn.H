{
    fvScalarMatrix hPorousEqn
    (
        fvm::div(porousPhi, porousH)
      - fvm::Sp(fvc::div(porousPhi), porousH)
      - fvm::laplacian(turbPorous.alphaEff(), porousH)
     ==
      - fvc::div(porousPhi, 0.5*magSqr(porousU), "div(phi,K)")
      + porousSources(porousRho, porousH)
    );

    hPorousEqn.relax();
    hPorousEqn.solve();

    porousThermo.correct();

    Info<< "Min/max in the porous T:"
        << min(porousThermo.T()).value() << ' '
        << max(porousThermo.T()).value() << endl;
}

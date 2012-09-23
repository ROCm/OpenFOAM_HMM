{
    volScalarField& hea = thermo.he();

    fvScalarMatrix EaEqn
    (
        fvm::ddt(rho, hea) + mvConvection->fvmDiv(phi, hea)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            hea.name() == "ea"
          ? fvc::div(phi, volScalarField("Ep", p/rho))
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), hea)
    );

    EaEqn.relax();
    EaEqn.solve();

    thermo.correct();
}

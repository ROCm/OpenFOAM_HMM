{
    volScalarField& hea = thermo.he();

    solve
    (
        betav*fvm::ddt(rho, hea) + mvConvection->fvmDiv(phi, hea)
      + betav*fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            hea.name() == "ea"
          ? fvc::div(phi, volScalarField("Ep", p/rho))
          : -betav*dpdt
        )
      - fvm::laplacian(Db, hea)
    );

    thermo.correct();
}

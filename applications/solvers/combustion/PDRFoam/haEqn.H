{
    solve
    (
        betav*fvm::ddt(rho, ha)
      + mvConvection->fvmDiv(phi, ha)
      - fvm::laplacian(Db, ha)
     ==
        betav*dpdt
      - betav*(fvc::ddt(rho, K) + fvc::div(phi, K))
    );

    thermo.correct();
}

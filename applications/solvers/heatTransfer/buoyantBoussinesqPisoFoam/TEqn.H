{
    volScalarField kappaEff
    (
        "kappaEff",
        turbulence->nu() + turbulence->nut()/Prt
    );

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
    );

    TEqn.relax();

    TEqn.solve();
}

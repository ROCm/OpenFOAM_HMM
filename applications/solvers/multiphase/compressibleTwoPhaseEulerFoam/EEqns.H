{
    volScalarField k1
    (
        "k1",
        alpha1*(thermo1.alpha()/rho1 + sqr(Ct)*nut2*thermo1.CpByCpv()/Prt)
    );

    volScalarField k2
    (
        "k2",
        alpha2*(thermo2.alpha()/rho2 + nut2*thermo2.CpByCpv()/Prt)
    );

    volScalarField& he1 = thermo1.he();
    volScalarField& he2 = thermo2.he();

    fvScalarMatrix he1Eqn
    (
        fvm::ddt(alpha1, he1)
      + fvm::div(alphaPhi1, he1)
      - fvm::laplacian(k1, he1)
     ==
        heatTransferCoeff*(thermo1.he(p, thermo2.T())/thermo1.Cp())/rho1
      - fvm::Sp(heatTransferCoeff/thermo1.Cp()/rho1, he1)
      + alpha1*(dpdt/rho1 - (fvc::ddt(K1) + fvc::div(phi1, K1)))
    );

    fvScalarMatrix he2Eqn
    (
        fvm::ddt(alpha2, he2)
      + fvm::div(alphaPhi2, he2)
      - fvm::laplacian(k2, he2)
     ==
        heatTransferCoeff*(thermo2.he(p, thermo1.T())/thermo2.Cp())/rho2
      - fvm::Sp(heatTransferCoeff/thermo2.Cp()/rho2, he2)
      + alpha2*(dpdt/rho2 - (fvc::ddt(K2) + fvc::div(phi2, K2)))
    );

    he1Eqn.relax();
    he1Eqn.solve();

    he2Eqn.relax();
    he2Eqn.solve();

    thermo1.correct();
    thermo2.correct();
}

if(turbulence)
{
    volScalarField divU = fvc::div(Uf & mesh.Sf());

    tmp<volTensorField> tgradU = fvc::grad(Uf);
    volScalarField G = 2*mut*(tgradU() && dev(symm(tgradU())));
    tgradU.clear();

    // Add the blockage generation term so that it is included consistently
    // in both the k and epsilon equations
    volScalarField GR = rho*mag(U)*(U & CT & U);

#   include "wallFunctions.H"

    // Dissipation equation
    fvScalarMatrix epsEqn
    (
        betav*fvm::ddt(rho, epsilon)
      + fvm::div(phi, epsilon)
      - fvm::laplacian(fvc::interpolate(alphaEps*muEff), epsilon)
      ==
        C1*(betav*G + GR)*epsilon/k
      - fvm::SuSp((2.0/3.0*C1)*betav*rho*divU, epsilon)
      - fvm::Sp(C2*betav*rho*epsilon/k, epsilon)
    );

#   include "wallDissipation.H"

    epsEqn.solve();
    bound(epsilon, dimensionedScalar("0", epsilon.dimensions(), 1.0e-15));


    // Turbulent kinetic energy equation
    solve
    (
        betav*fvm::ddt(rho, k)
      + fvm::div(phi, k)
      - fvm::laplacian(fvc::interpolate(alphak*muEff), k)
      ==
        betav*G + GR
      - fvm::SuSp(2.0/3.0*betav*rho*divU, k)
      - fvm::Sp(betav*rho*epsilon/k, k)
    );

    bound(k, dimensionedScalar("0", k.dimensions(), 0.0));

    //- Re-calculate turbulence viscosity
    mut = Cmu*rho*sqr(k)/epsilon;

#   include "wallViscosity.H"

}

muEff = mut + thermo->mu();

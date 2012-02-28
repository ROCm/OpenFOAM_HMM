if (turbulence)
{
    if (mesh.changing())
    {
        y.correct();
    }

    tmp<volTensorField> tgradUb = fvc::grad(Ub);
    volScalarField G(2*nutb*(tgradUb() && dev(symm(tgradUb()))));
    tgradUb.clear();

    #include "wallFunctions.H"

    // Dissipation equation
    fvScalarMatrix epsEqn
    (
        fvm::ddt(epsilon)
      + fvm::div(phib, epsilon)
      - fvm::Sp(fvc::div(phib), epsilon)
      - fvm::laplacian
        (
            alphaEps*nuEffb, epsilon,
            "laplacian(DepsilonEff,epsilon)"
        )
      ==
         C1*G*epsilon/k
       - fvm::Sp(C2*epsilon/k, epsilon)
    );

    #include "wallDissipation.H"

    epsEqn.relax();
    epsEqn.solve();

    epsilon.max(dimensionedScalar("zero", epsilon.dimensions(), 1.0e-15));


    // Turbulent kinetic energy equation
    fvScalarMatrix kEqn
    (
        fvm::ddt(k)
      + fvm::div(phib, k)
      - fvm::Sp(fvc::div(phib), k)
      - fvm::laplacian
        (
            alphak*nuEffb, k,
            "laplacian(DkEff,k)"
        )
      ==
        G
      - fvm::Sp(epsilon/k, k)
    );
    kEqn.relax();
    kEqn.solve();

    k.max(dimensionedScalar("zero", k.dimensions(), 1.0e-8));

    //- Re-calculate turbulence viscosity
    nutb = Cmu*sqr(k)/epsilon;

    #include "wallViscosity.H"
}

nuEffb = nutb + nub;

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho.dimensionedInternalField()*g
      + reactingParcels.SU()
    );

    UEqn.relax();

    tmp<volScalarField> trAU;
    tmp<volTensorField> trTU;

    if (pressureImplicitPorosity)
    {
        tmp<volTensorField> tTU = tensor(I)*UEqn.A();
        pZones.addResistance(UEqn, tTU());
        trTU = inv(tTU());
        trTU().rename("rAU");

        volVectorField gradp = fvc::grad(p);

        for (int UCorr=0; UCorr<nUCorr; UCorr++)
        {
            U = trTU() & (UEqn.H() - gradp);
        }
        U.correctBoundaryConditions();
    }
    else
    {
        pZones.addResistance(UEqn);

        solve
        (
            UEqn == -fvc::grad(p)
        );

        trAU = 1.0/UEqn.A();
        trAU().rename("rAU");
    }

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }

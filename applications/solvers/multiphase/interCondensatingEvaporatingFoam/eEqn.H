{
    tmp<volScalarField> tcp(thermo->Cp());
    const volScalarField& cp = tcp();

    kappaEff = thermo->kappa() + rho*cp*turbulence->nut()/Prt;

    pDivU = dimensionedScalar("pDivU", p.dimensions()/dimTime, 0.0);

    if (thermo->pDivU())
    {
        pDivU = (p*fvc::div(rhoPhi/fvc::interpolate(rho)));
    }

    fvScalarMatrix eEqn
    (
        fvm::ddt(rho, e)
      + fvm::div(rhoPhi, e)
      - fvm::laplacian(kappaEff/cp, e)
      + pDivU
    );

    eEqn.relax();
    eEqn.solve();

    thermo->correct();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}

{
    tmp<volScalarField> pWork
    (
        new volScalarField
        (
            IOobject
            (
                "pWork",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );

    if (pressureWork)
    {
        surfaceScalarField phiU("phiU", phi/fvc::interpolate(rho));

        pWork() += fvc::div(phiU*fvc::interpolate(p)) - p*fvc::div(phiU);

        if (pressureWorkTimeDerivative)
        {
            pWork() += fvc::ddt(p);
        }
    }

    {
        solve
        (
            fvm::ddt(rho, hs)
          + mvConvection->fvmDiv(phi, hs)
          - fvm::laplacian(turbulence->alphaEff(), hs)
         ==
            pWork()
          + parcels.Sh(hs)
          + radiation->Shs(thermo)
          + energySource.Su()
          + chemistrySh
        );

        thermo.correct();

        radiation->correct();

        Info<< "T gas min/max   = " << min(T).value() << ", "
            << max(T).value() << endl;
    }
}

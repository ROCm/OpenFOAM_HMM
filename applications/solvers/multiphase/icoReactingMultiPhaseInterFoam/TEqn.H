{
    radiation->correct();

    tmp<volScalarField> texpSource
    (
        new volScalarField
        (
            IOobject
            (
                "texpSource",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimTemperature/dimTime, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& expSource = texpSource.ref();

    tmp<volScalarField> tkappaEff
    (
        new volScalarField
        (
            IOobject
            (
                "kappaEff",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", sqr(dimLength)/dimTime, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& kappaEff = tkappaEff.ref();

    //const surfaceScalarField rhoTempPhi("phi", rhoPhi/fvc::interpolate(rho));

    const surfaceScalarField rhoTempPhi("phi", fluid.phi());

    const volScalarField divU(fvc::div(phi));

    forAllIter(UPtrList<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();
        const volScalarField& alpha = phase;
        const volScalarField DDtAlpha(fvc::DDt(phi, alpha));

        const volScalarField invCpRho(1.0/phase.rho()/phase.Cp());

        if (fluid.dpdt())
        {
            const volScalarField ddtp(fvc::ddt(p));
            expSource += (DDtAlpha*p + alpha*(p*divU + ddtp))*invCpRho;
        }
        else
        {
            expSource += (DDtAlpha*p + alpha*(p*divU))*invCpRho;
        }

        kappaEff += alpha*phase.kappa()*invCpRho;

        DebugVar(max(alpha*phase.kappa()));
        DebugVar(max(alpha*phase.Cp()));
    }

    kappaEff += turbulence->nut()/fluid.Prt();

    if (mesh.time().outputTime())
    {
        expSource.write();
        kappaEff.write();
    }

    //dimensionedScalar S("S", dimEnergy/dimVolume/dimTime, 1.225e8);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(rhoTempPhi, T)
      - fvm::laplacian(kappaEff, T,  "laplacian(kappa,T)")
      ==
        fluid.heatTransfer(T)
      + expSource
      + radiation->ST(fluid.Cp()*rho, T)
//      + S/Cp/rho
    );

    TEqn.relax();
    TEqn.solve();

    fluid.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}

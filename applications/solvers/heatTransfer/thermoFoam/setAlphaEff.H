    Info<< "Creating turbulence model\n" << endl;
    tmp<volScalarField> talphaEff;

    IOobject turbulenceHeader
    (
        "turbulenceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject RASHeader
    (
        "RASProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject LESHeader
    (
        "LESProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ
    );

    if (turbulenceHeader.typeHeaderOk<IOdictionary>(true))
    {
        autoPtr<compressible::turbulenceModel> turbulence
        (
            compressible::turbulenceModel::New
            (
                rho,
                U,
                phi,
                thermo
            )
        );

        talphaEff = turbulence->alphaEff();
    }
    else if (RASHeader.typeHeaderOk<IOdictionary>(true))
    {
        autoPtr<compressible::RASModel> turbulence
        (
            compressible::RASModel::New
            (
                rho,
                U,
                phi,
                thermo
            )
        );

        talphaEff = turbulence->alphaEff();
    }
    else if (LESHeader.typeHeaderOk<IOdictionary>(true))
    {
        autoPtr<compressible::LESModel> turbulence
        (
            compressible::LESModel::New
            (
                rho,
                U,
                phi,
                thermo
            )
        );

        talphaEff = turbulence->alphaEff();
    }
    else
    {
        talphaEff = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "alphaEff",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimMass/dimLength/dimTime, 0.0)
            )
        );
    }

    const volScalarField& alphaEff = talphaEff();

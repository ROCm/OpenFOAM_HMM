    autoPtr<volScalarField> allRho;
    autoPtr<volScalarField> allh;
    autoPtr<volScalarField> allVolAlpha;
    autoPtr<fvScalarMatrix> allSource;

    if (temperatureCoupled)
    {
        allRho.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    allMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                allMesh(),
                dimensionedScalar("rho", dimDensity, 0.0)
            )
        );

        allh.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "h",
                    runTime.timeName(),
                    allMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                allMesh(),
                dimensionedScalar("h", dimEnergy/dimMass, 0.0),
                mixedFvPatchScalarField::typeName
            )
        );

        allVolAlpha.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "volAlpha",
                    runTime.timeName(),
                    allMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                allMesh(),
                dimensionedScalar("volAlpha", dimMass/dimLength/dimTime, 0.0)
            )
        );

        allSource.reset
        (
            new fvMatrix<scalar>(allh(), allh().dimensions()*dimMass/dimTime)
        );
    }

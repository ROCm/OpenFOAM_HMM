    label nAveragingSteps = 1;

    IFstream nAveragingStepsFile
    (
        runTime.path()/runTime.timeName()/"uniform"/"nAveragingSteps.raw"
    );


    autoPtr<volVectorField> UmeanPtr;
    autoPtr<volSymmTensorField> RPtr;
    autoPtr<volSymmTensorField> BmeanPtr;
    autoPtr<volScalarField> epsilonMeanPtr;

    autoPtr<volScalarField> pMeanPtr;
    autoPtr<volScalarField> pPrime2MeanPtr;

    if (nAveragingStepsFile.good())
    {
        nAveragingStepsFile >> nAveragingSteps;

        Info<< "Reading field Umean\n" << endl;
        UmeanPtr.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "Umean",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        Info<< "\nReading field R\n" << endl;
        RPtr.reset
        (
            new volSymmTensorField
            (
                IOobject
                (
                    "R",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        Info<< "\nReading field Bmean\n" << endl;
        BmeanPtr.reset
        (
            new volSymmTensorField
            (
                IOobject
                (
                    "Bmean",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        Info<< "\nReading field epsilonMean\n" << endl;
        epsilonMeanPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "epsilonMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        Info<< "Reading field pMean\n" << endl;
        pMeanPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "pMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        Info<< "Reading field pPrime2Mean\n" << endl;
        pPrime2MeanPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "pPrime2Mean",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "Creating field Umean\n" << endl;
        UmeanPtr.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "Umean",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                U
            )
        );

        Info<< "Creating field R\n" << endl;
        RPtr.reset
        (
            new volSymmTensorField
            (
                IOobject
                (
                    "R",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                sqr(U) - sqr(UmeanPtr())
            )
        );

        Info<< "Creating field Bmean\n" << endl;
        BmeanPtr.reset
        (
            new volSymmTensorField
            (
                IOobject
                (
                    "Bmean",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                sgsModel->B()
            )
        );

        Info<< "Creating field epsilonMean\n" << endl;
        epsilonMeanPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "epsilonMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                sgsModel->epsilon()
            )
        );

        Info<< "Creating field pMean\n" << endl;
        pMeanPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "pMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                p
            )
        );

        Info<< "Creating field pPrime2Mean\n" << endl;
        pPrime2MeanPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "pPrime2Mean",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                sqr(p) - sqr(pMeanPtr())
            )
        );
    }

    volVectorField& Umean = UmeanPtr();
    volSymmTensorField& R = RPtr();
    volSymmTensorField& Bmean = BmeanPtr();
    volScalarField& epsilonMean = epsilonMeanPtr();

    volScalarField& pMean = pMeanPtr();
    volScalarField& pPrime2Mean = pPrime2MeanPtr();

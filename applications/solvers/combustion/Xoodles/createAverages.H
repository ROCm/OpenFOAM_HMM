    label nAveragingSteps = 1;

    IFstream nAveragingStepsFile
    (
        runTime.path()/runTime.timeName()/"uniform"/"nAveragingSteps.raw"
    );

    volVectorField* UmeanPtr;
    volSymmTensorField* RPtr;

    volScalarField* pMeanPtr;
    volScalarField* pPrime2MeanPtr;

    volScalarField* TmeanPtr;
    volScalarField* bmeanPtr;
    volScalarField* XimeanPtr;

    if (nAveragingStepsFile.good())
    {
        nAveragingStepsFile >> nAveragingSteps;

        Info << "Reading field Umean\n" << endl;
        UmeanPtr = new volVectorField
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
        );

        Info<< "\nReading field R\n" << endl;
        RPtr = new volSymmTensorField
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
        );

        Info<< "Reading field pMean\n" << endl;
        pMeanPtr = new volScalarField
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
        );

        Info<< "Reading field pPrime2Mean\n" << endl;
        pPrime2MeanPtr = new volScalarField
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
        );

        Info << "Reading field Tmean\n" << endl;
        TmeanPtr = new volScalarField
        (
            IOobject
            (
                "Tmean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info << "Reading field bmean\n" << endl;
        bmeanPtr = new volScalarField
        (
            IOobject
            (
                "bmean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info << "Reading field Ximean\n" << endl;
        XimeanPtr = new volScalarField
        (
            IOobject
            (
                "Ximean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
    }
    else
    {
        Info << "Creating field Umean\n" << endl;
        UmeanPtr = new volVectorField
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
        );

        Info<< "Creating field R\n" << endl;
        RPtr = new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqr(U) - sqr(*UmeanPtr)
        );

        Info<< "Creating field pMean\n" << endl;
        pMeanPtr = new volScalarField
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
        );

        Info<< "Creating field pPrime2Mean\n" << endl;
        pPrime2MeanPtr = new volScalarField
        (
            IOobject
            (
                "pPrime2Mean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqr(p) - sqr(*pMeanPtr)
        );

        Info << "Creating field Tmean\n" << endl;
        TmeanPtr = new volScalarField
        (
            IOobject
            (
                "Tmean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            T
        );

        Info << "Creating field bmean\n" << endl;
        bmeanPtr = new volScalarField
        (
            IOobject
            (
                "bmean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            b
        );

        Info << "Creating field Ximean\n" << endl;
        XimeanPtr = new volScalarField
        (
            IOobject
            (
                "Ximean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Xi
        );
    }

    volVectorField& Umean = *UmeanPtr;
    volSymmTensorField& R = *RPtr;

    volScalarField& pMean = *pMeanPtr;
    volScalarField& pPrime2Mean = *pPrime2MeanPtr;

    volScalarField& Tmean = *TmeanPtr;
    volScalarField& bmean = *bmeanPtr;
    volScalarField& Ximean = *XimeanPtr;

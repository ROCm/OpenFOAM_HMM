    label nAveragingSteps = 1;

    IFstream nAveragingStepsFile
    (
        runTime.path()/runTime.timeName()/"uniform"/"nAveragingSteps.raw"
    );


    volVectorField* UmeanPtr;
    //volTensorField* RPtr;
    //volTensorField* BmeanPtr;
    volScalarField* kMeanPtr;
    volScalarField* kGSMeanPtr;
    volScalarField* nuEffMeanPtr;
    volScalarField* gammaMeanPtr;
    volScalarField* sigmaMeanPtr;

    volScalarField* pMeanPtr;
    volScalarField* pPrime2MeanPtr;

    if (nAveragingStepsFile.good())
    {
        nAveragingStepsFile >> nAveragingSteps;

        Info<< "Reading field Umean\n" << endl;
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
        /*
        Info<< "\nReading field R\n" << endl;
        RPtr = new volTensorField
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

        Info<< "\nReading field Bmean\n" << endl;
        BmeanPtr = new volTensorField
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
        );
        */
        Info<< "\nReading field kMean\n" << endl;
        kMeanPtr = new volScalarField
        (
            IOobject
            (
                "kMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
        Info<< "\nReading field kGSMean\n" << endl;
        kGSMeanPtr = new volScalarField
        (
            IOobject
            (
                "kGSMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );


        Info<< "\nReading field nuEffMean\n" << endl;
        nuEffMeanPtr = new volScalarField
        (
            IOobject
            (
                "nuEffMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
        Info<< "\nReading field gammaMean\n" << endl;
        gammaMeanPtr = new volScalarField
        (
            IOobject
            (
                "gammaMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
        Info<< "\nReading field sigmaMean\n" << endl;
        sigmaMeanPtr = new volScalarField
        (
            IOobject
            (
                "sigmaMean",
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
    }
    else
    {
        Info<< "Creating field Umean\n" << endl;
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
        /*
        Info<< "Creating field R\n" << endl;
        RPtr = new volTensorField
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

        Info<< "Creating field Bmean\n" << endl;
        BmeanPtr = new volTensorField
        (
            IOobject
            (
                "Bmean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            turbulence->B()
        );
        */
        Info<< "Creating field kMean\n" << endl;
        kMeanPtr = new volScalarField
        (
            IOobject
            (
                "kMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            turbulence->k()
        );
        Info<< "Creating field kGSMean\n" << endl;
        kGSMeanPtr = new volScalarField
        (
            IOobject
            (
                "kGSMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            0.5*tr(sqr(U) - sqr(*UmeanPtr))
        );
    
    
    
        Info<< "Creating field nuEffMean\n" << endl;
        nuEffMeanPtr = new volScalarField
        (
            IOobject
            (
                "nuEffMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            turbulence->nuEff()
        );
        Info<< "Creating field gammaMean\n" << endl;
        gammaMeanPtr = new volScalarField
        (
            IOobject
            (
                "gammaMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gamma
        );
        Info<< "Creating field sigmaMean\n" << endl;
        sigmaMeanPtr = new volScalarField
        (
            IOobject
            (
                "sigmaMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mag(fvc::grad(gamma))
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
            pd
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
            sqr(pd) - sqr(*pMeanPtr)
        );
    }

    volVectorField& Umean = *UmeanPtr;
    //volTensorField& R = *RPtr;
    //volTensorField& Bmean = *BmeanPtr;
    volScalarField& kMean = *kMeanPtr;
    volScalarField& kGSMean = *kGSMeanPtr;
    volScalarField& nuEffMean = *nuEffMeanPtr;
    volScalarField& gammaMean = *gammaMeanPtr;
    volScalarField& sigmaMean = *sigmaMeanPtr;

    volScalarField& pMean = *pMeanPtr;
    volScalarField& pPrime2Mean = *pPrime2MeanPtr;

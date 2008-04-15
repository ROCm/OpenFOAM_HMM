    label nAveragingSteps = 1;

    IFstream nAveragingStepsFile
    (
        runTime.path()/runTime.timeName()/"nAveragingSteps.raw"
    );


    volVectorField* UaMeanPtr;
    volVectorField* UbMeanPtr;
 
    volScalarField* alphaMeanPtr;
    volScalarField* pMeanPtr;

    if (nAveragingStepsFile.good())
    {
        nAveragingStepsFile >> nAveragingSteps;

        Info<< "Reading field UaMean\n" << endl;
        UaMeanPtr = new volVectorField
        (
            IOobject
            (
                "UaMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
        
        Info<< "Reading field UbMean\n" << endl;
        UbMeanPtr = new volVectorField
        (
            IOobject
            (
                "UbMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info<< "Reading field alphaMean\n" << endl;
        alphaMeanPtr = new volScalarField
        (
            IOobject
            (
                "alphaMean",
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
    }
    else
    {
        Info<< "Creating field UaMean\n" << endl;
        UaMeanPtr = new volVectorField
        (
            IOobject
            (
                "UaMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Ua
        );
        
        Info<< "Creating field UbMean\n" << endl;
        UbMeanPtr = new volVectorField
        (
            IOobject
            (
                "UbMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Ub
        );

        Info<< "Creating field alphaMean\n" << endl;
        alphaMeanPtr = new volScalarField
        (
            IOobject
            (
                "alphaMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            alpha
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
    }

    volVectorField& UaMean = *UaMeanPtr;
    volVectorField& UbMean = *UbMeanPtr;

    volScalarField& alphaMean = *alphaMeanPtr;
    volScalarField& pMean = *pMeanPtr;

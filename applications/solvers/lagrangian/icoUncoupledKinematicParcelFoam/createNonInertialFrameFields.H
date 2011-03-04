    Info<< "Reading non-inertial frame fields" << endl;

    IOobject linearAccelerationHeader
    (
        "linearAcceleration",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> linearAccelerationPtr;

    if (linearAccelerationHeader.headerOk())
    {
        Info<< "    Reading " << linearAccelerationHeader.name() << endl;

        linearAccelerationPtr.reset
        (
            new uniformDimensionedVectorField(linearAccelerationHeader)
        );
    }


    IOobject angularVelocityMagHeader
    (
        "angularVelocityMag",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedScalarField> angularVelocityMagPtr;

    if (angularVelocityMagHeader.headerOk())
    {
        Info<< "    Reading " << angularVelocityMagHeader.name() << endl;

        angularVelocityMagPtr.reset
        (
            new uniformDimensionedScalarField(angularVelocityMagHeader)
        );
    }


    IOobject angularAccelerationHeader
    (
        "angularAcceleration",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> angularAccelerationPtr;

    if (angularAccelerationHeader.headerOk())
    {
        Info<< "    Reading " << angularAccelerationHeader.name() << endl;

        angularAccelerationPtr.reset
        (
            new uniformDimensionedVectorField(angularAccelerationHeader)
        );
    }


    IOobject axisHeader
    (
        "axis",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> axisPtr;

    if (axisHeader.headerOk())
    {
        Info<< "    Reading " << axisHeader.name() << endl;

        axisPtr.reset
        (
            new uniformDimensionedVectorField(axisHeader)
        );
    }


    IOobject axisRefPointHeader
    (
        "axisRefPoint",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    autoPtr<uniformDimensionedVectorField> axisRefPointPtr;

    if (axisRefPointHeader.headerOk())
    {
        Info<< "    Reading " << axisRefPointHeader.name() << endl;

        axisRefPointPtr.reset
        (
            new uniformDimensionedVectorField(axisRefPointHeader)
        );
    }

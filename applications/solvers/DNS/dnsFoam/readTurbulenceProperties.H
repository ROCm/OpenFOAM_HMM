    Info<< "Reading turbulenceProperties\n" << endl;

    IOdictionary turbulenceProperties
    (
        IOobject
        (
            "turbulenceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    volVectorField force
    (
        U/dimensionedScalar("dt", dimTime, runTime.deltaTValue())
    );

    Kmesh K(mesh);
    UOprocess forceGen(K, runTime.deltaTValue(), turbulenceProperties);

    label ntot = 1;
    forAll(K.nn(), idim)
    {
        ntot *= K.nn()[idim];
    }
    const scalar recRootN = 1.0/Foam::sqrt(scalar(ntot));

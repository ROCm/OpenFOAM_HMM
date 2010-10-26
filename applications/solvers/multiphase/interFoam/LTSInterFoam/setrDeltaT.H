{
    scalar maxCo
    (
        piso.lookupOrDefault<scalar>("maxCo", 0.9)
    );

    scalar maxAlphaCo
    (
        piso.lookupOrDefault<scalar>("maxAlphaCo", 0.2)
    );

    scalar rDeltaTSmoothingCoeff
    (
        piso.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    label nAlphaSpreadIter
    (
        piso.lookupOrDefault<label>("nAlphaSpreadIter", 1)
    );

    label nAlphaSweepIter
    (
        piso.lookupOrDefault<label>("nAlphaSweepIter", 5)
    );

    scalar rDeltaTDampingCoeff
    (
        piso.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    scalar maxDeltaT
    (
        piso.lookupOrDefault<scalar>("maxDeltaT", GREAT)
    );

    volScalarField rDeltaT0 = rDeltaT;

    // Set the reciprocal time-step using an effective maximum Courant number
    rDeltaT = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum
        (
            mag(rhoPhi)*mesh.deltaCoeffs()/(maxCo*mesh.magSf())
        )/rho
    );

    // Limit the time-step further in the region of the interface
    {
        surfaceScalarField alphaf = fvc::interpolate(alpha1);

        surfaceScalarField SfUfbyDelta =
            pos(alphaf - 0.01)*pos(0.99 - alphaf)
           *mesh.surfaceInterpolation::deltaCoeffs()*mag(phi);

        rDeltaT = max
        (
            rDeltaT,
            fvc::surfaceSum(mag(SfUfbyDelta/(maxAlphaCo*mesh.magSf())))
        );
    }
    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.internalField())
        << ", " << gMax(1/rDeltaT.internalField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    if (nAlphaSpreadIter > 0)
    {
        fvc::spread(rDeltaT, alpha1, nAlphaSpreadIter);
    }

    if (nAlphaSweepIter > 0)
    {
        fvc::sweep(rDeltaT, alpha1, nAlphaSweepIter);
    }

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.internalField())
        << ", " << gMax(1/rDeltaT.internalField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        Info<< "Damping rDeltaT" << endl;
        rDeltaT = rDeltaT0*max(rDeltaT/rDeltaT0, 1.0 - rDeltaTDampingCoeff);
    }

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.internalField())
        << ", " << gMax(1/rDeltaT.internalField()) << endl;

    label nAlphaSubCycles
    (
        readLabel(piso.lookup("nAlphaSubCycles"))
    );

    rSubDeltaT = rDeltaT*nAlphaSubCycles;
}

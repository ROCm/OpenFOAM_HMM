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

    scalar alphaSpreadDiff
    (
        piso.lookupOrDefault<scalar>("alphaSpreadDiff", 0.2)
    );

    scalar alphaSpreadMax
    (
        piso.lookupOrDefault<scalar>("alphaSpreadMax", 0.99)
    );

    scalar alphaSpreadMin
    (
        piso.lookupOrDefault<scalar>("alphaSpreadMin", 0.01)
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

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.dimensionedInternalField() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(mag(rhoPhi))().dimensionedInternalField()
       /((2*maxCo)*mesh.V()*rho.dimensionedInternalField())
    );

    if (maxAlphaCo < maxCo)
    {
        // Further limit the reciprocal time-step
        // in the vicinity of the interface

        volScalarField alpha1Bar = fvc::average(alpha1);

        rDeltaT.dimensionedInternalField() = max
        (
            rDeltaT.dimensionedInternalField(),
            pos(alpha1Bar.dimensionedInternalField() - alphaSpreadMin)
           *pos(alphaSpreadMax - alpha1Bar.dimensionedInternalField())
           *fvc::surfaceSum(mag(phi))().dimensionedInternalField()
           /((2*maxAlphaCo)*mesh.V())
        );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.internalField())
        << ", " << gMax(1/rDeltaT.internalField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    if (nAlphaSpreadIter > 0)
    {
        fvc::spread
        (
            rDeltaT,
            alpha1,
            nAlphaSpreadIter,
            alphaSpreadDiff,
            alphaSpreadMax,
            alphaSpreadMin
        );
    }

    if (nAlphaSweepIter > 0)
    {
        fvc::sweep(rDeltaT, alpha1, nAlphaSweepIter, alphaSpreadDiff);
    }

    Info<< "Smoothed flow time scale min/max = "
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
        rDeltaT =
             rDeltaT0
            *max(rDeltaT/rDeltaT0, scalar(1.0) - rDeltaTDampingCoeff);

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.internalField())
            << ", " << gMax(1/rDeltaT.internalField()) << endl;
    }

    label nAlphaSubCycles
    (
        readLabel(piso.lookup("nAlphaSubCycles"))
    );

    rSubDeltaT = rDeltaT*nAlphaSubCycles;
}

    scalar CoNum = -GREAT;
    forAll(fluidRegions, regionI)
    {
        CoNum = Foam::max
        (
            compressibleCourantNo
            (
                fluidRegions[regionI],
                runTime,
                rhoFluid[regionI],
                phiFluid[regionI]
            ),
            CoNum
        );
    }

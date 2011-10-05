    scalar CoNum = -GREAT;

    forAll(fluidToProc, regionI)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                procMeshes[fluidToProc[regionI]],
                runTime,
                rhoFluid[regionI],
                phiFluid[regionI]
            ),
            CoNum
        );
    }

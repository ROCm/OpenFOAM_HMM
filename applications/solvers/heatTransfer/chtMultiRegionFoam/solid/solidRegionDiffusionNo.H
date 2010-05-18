    scalar DiNum = -GREAT;

    forAll(solidRegions, regionI)
    {
        DiNum = max
        (
            solidRegionDiffNo
            (
                solidRegions[regionI],
                runTime,
                rhos[regionI]*cps[regionI],
                Ks[regionI]
            ),
            DiNum
        );
    }

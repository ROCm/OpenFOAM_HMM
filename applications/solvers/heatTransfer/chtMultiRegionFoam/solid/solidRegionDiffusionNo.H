    scalar DiNum = -GREAT;
    if (solidRegions.size())//
    {
        forAll(solidRegions, regionI)
        {
            DiNum = max
            (
                solidRegionDiffNo
                (
                    solidRegions[regionI],
                    runTime,
                    rhosCps[regionI],
                    Ks[regionI]
                ),
                DiNum
            );
        }
    }
